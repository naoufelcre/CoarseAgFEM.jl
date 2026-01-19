module RobustAgFEM

using Gridap
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Helpers
using GridapEmbedded
using GridapEmbedded.LevelSetCutters
using GridapEmbedded.CSG
using GridapEmbedded.Interfaces: CUT, IN, OUT, CUT_IN, CUT_OUT
using LinearAlgebra

export RobustAggregation

"""
    RobustAggregation(threshold=0.5, max_iters=100)

Strategy for AgFEM aggregation that allows configuring the maximum number of iterations
to handle long chains of small cut elements often found in adaptive meshes.
"""
struct RobustAggregation
    threshold::Float64
    max_iters::Int
end

RobustAggregation(t) = RobustAggregation(t, 200)

function GridapEmbedded.AgFEM.aggregate(
    strategy::RobustAggregation,
    cut::EmbeddedDiscretization,
    geo::CSG.Geometry,
    in_or_out)

    facet_to_inoutcut = GridapEmbedded.AgFEM.compute_bgfacet_to_inoutcut(cut.bgmodel,geo)
    _robust_aggregate_by_threshold(
        strategy.threshold,
        strategy.max_iters,
        cut,geo,in_or_out,facet_to_inoutcut
    )
end

function _robust_aggregate_by_threshold(threshold,max_iters,cut,geo,loc,facet_to_inoutcut)
    @assert loc in (IN,OUT)

    cutinorout = loc == IN ? (CUT_IN,IN) : (CUT_OUT,OUT)
    trian = Gridap.Geometry.Triangulation(cut,cutinorout,geo)
    model = cut.bgmodel
    bgtrian = Gridap.Geometry.get_triangulation(model)
    cell_to_cut_meas = Gridap.CellData.get_cell_measure(trian,bgtrian)
    cell_to_meas = Gridap.CellData.get_cell_measure(bgtrian)
    cell_to_unit_cut_meas = Gridap.Arrays.lazy_map(/,cell_to_cut_meas,cell_to_meas)

    cell_to_inoutcut = GridapEmbedded.AgFEM.compute_bgcell_to_inoutcut(cut,geo)

    cell_to_coords = Gridap.Geometry.get_cell_coordinates(bgtrian)
    topo = Gridap.Geometry.get_grid_topology(model)
    D = Gridap.Geometry.num_cell_dims(model)
    cell_to_faces = Gridap.Geometry.get_faces(topo,D,D-1)
    face_to_cells = Gridap.Geometry.get_faces(topo,D-1,D)

    _robust_aggregate_kernel(
        threshold, max_iters,
        cell_to_unit_cut_meas,facet_to_inoutcut,cell_to_inoutcut,
        loc,cell_to_coords,cell_to_faces,face_to_cells)
end

function _robust_aggregate_kernel(
    threshold, max_iters,
    cell_to_unit_cut_meas,facet_to_inoutcut,cell_to_inoutcut,
    loc,cell_to_coords,cell_to_faces,face_to_cells)

    n_cells = length(cell_to_unit_cut_meas)
    cell_to_cellin = zeros(Int32,n_cells)
    cell_to_touched = fill(false,n_cells)

    # 1. Initialize Stable Cells
    for cell in 1:n_cells
        if cell_to_unit_cut_meas[cell] ≥ threshold
            cell_to_cellin[cell] = cell
            cell_to_touched[cell] = true
        end
    end

    c1 = Gridap.Arrays.array_cache(cell_to_faces)
    c2 = Gridap.Arrays.array_cache(face_to_cells)
    c3 = Gridap.Arrays.array_cache(cell_to_coords)
    c4 = Gridap.Arrays.array_cache(cell_to_coords)

    all_aggregated = false
    
    # 2. Iterative Aggregation
    for iter in 1:max_iters
        all_aggregated = true
        for cell in 1:n_cells
            # If cell needs aggregation (not touched, and is CUT)
            if ! cell_to_touched[cell] && cell_to_inoutcut[cell] == CUT
                neigh_cell = _robust_find_best_neighbor(
                    c1,c2,c3,c4,cell,
                    cell_to_faces,
                    face_to_cells,
                    cell_to_coords,
                    cell_to_touched,
                    cell_to_cellin,
                    facet_to_inoutcut,
                    loc)
                
                if neigh_cell > 0
                    cellin = cell_to_cellin[neigh_cell]
                    cell_to_cellin[cell] = cellin
                else
                    all_aggregated = false
                end
            end
        end
        if all_aggregated
            println("    [RobustAgFEM] Aggregated all cells in $iter iterations.")
            break
        end
        
        # Update touched status for next pass
        _robust_touch_aggregated_cells!(cell_to_touched,cell_to_cellin)
    end
    

    if !all_aggregated
        warn_msg = "[RobustAgFEM] Failed to aggregate all cells after $max_iters iterations."
        println(warn_msg)
        # We do NOT error, but we return the partial aggregation. 
        # GridapEmbedded usually asserts. We'll warn.
        # Check how many failed
        n_failed = count(c -> (!cell_to_touched[c] && cell_to_inoutcut[c] == CUT), 1:n_cells)
        println("              -> $n_failed cells orphaned.")
    end

    cell_to_cellin
end

function _robust_find_best_neighbor(
    c1,c2,c3,c4,cell,
    cell_to_faces,
    face_to_cells,
    cell_to_coords,
    cell_to_touched,
    cell_to_cellin,
    facet_to_inoutcut,
    loc)

    faces = Gridap.Arrays.getindex!(c1,cell_to_faces,cell)
    dmin = Inf
    T = eltype(eltype(face_to_cells))
    best_neigh_cell = zero(T)
    i_to_coords = Gridap.Arrays.getindex!(c3,cell_to_coords,cell)
    
    for face in faces
        inoutcut = facet_to_inoutcut[face]
        # Only cross faces that are CUT or inside the fluid domain (LOC)
        if inoutcut != CUT && inoutcut != loc
            continue
        end
        
        neighs = Gridap.Arrays.getindex!(c2,face_to_cells,face)
        for neigh_cell in neighs
            # Determine if neighbor is a valid root candidate
            # Must be different cell, and already part of an aggregate (touched)
            if neigh_cell != cell && cell_to_touched[neigh_cell]
                cellin = cell_to_cellin[neigh_cell]
                j_to_coords = Gridap.Arrays.getindex!(c4,cell_to_coords,cellin)
                
                # Metric: Max distance between nodes of current cell and root cell
                # This penalizes long chains spatially
                d = 0.0
                for p in i_to_coords
                    for q in j_to_coords
                        d = max(d,Float64(norm(p-q)))
                    end
                end
                
                if (1.0+1.0e-9)*d < dmin
                    dmin = d
                    best_neigh_cell = neigh_cell
                end
            end
        end
    end
    best_neigh_cell
end

function _robust_touch_aggregated_cells!(cell_to_touched,cell_to_cellin)
    for (cell,cellin) in enumerate(cell_to_cellin)
        if cellin > 0
            cell_to_touched[cell] = true
        end
    end
end

end # module
