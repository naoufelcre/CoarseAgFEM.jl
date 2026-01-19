module GridapIntegration

using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using ..QuadtreeAggregations

export quadtree_to_discrete_model

"""
    quadtree_to_discrete_model(elements::Vector{QuadElement})

Converts a list of `QuadElement` (from QuadtreeAggregations) into a `Gridap.UnstructuredDiscreteModel`.
Handles mixed meshes (Triangles and Quadrilaterals).
"""
function quadtree_to_discrete_model(elements::Vector{QuadElement})
    println("  [GridapBridge] Converting $(length(elements)) elements to Gridap Model...")

    # 1. Extract Unique Nodes
    # -----------------------
    # Robust Merging: Tolerance based
    # Key is Tuple for hashing robustness
    # We use a coarse rounding for the Dict to bucket candidates, but check distance
    
    # Robust Merging: Rounding-based bucket
    # Key is Tuple (x, y) rounded
    node_to_id = Dict{Tuple{Float64, Float64}, Int}()
    node_coords = Vector{Point{2,Float64}}()
    
    n_queries = 0
    n_merged = 0
    
    function get_node_id(p::Vector{Float64})
        n_queries += 1
        # Round to 10 digits (~1e-10 precision)
        # This acts as a robust grid snap
        k = (round(p[1], digits=10), round(p[2], digits=10))
        
        if haskey(node_to_id, k)
            n_merged += 1
            return node_to_id[k]
        end
        
        push!(node_coords, Point(p[1], p[2]))
        id = length(node_coords)
        node_to_id[k] = id
        return id
    end

    # 2. Build Cell Connectivity & Types
    # ----------------------------------
    # We support MIXED meshes (TRI and QUAD)
    cell_node_ids = Vector{Vector{Int}}()
    cell_types = Int8[] 
    
    # Map to track lineage: Leaf ID -> [Cell ID 1, Cell ID 2, ...]
    leaf_to_cell_ids = Dict{Int, Vector{Int}}()

    # Counter for current cell index
    current_cell_idx = 0

    # Sort elements by type (Triangles first, then Quads) to help Gridap VTK writer
    # which might struggle with interleaved mixed strides
    sort!(elements, by = x -> length(x.nodes))

    for el in elements
        ids = [get_node_id(p) for p in el.nodes]
        
        created_cell_indices = Int[]

        if length(ids) == 3
            push!(cell_node_ids, ids)
            push!(cell_types, 1) # 1 = TRI
            current_cell_idx += 1
            push!(created_cell_indices, current_cell_idx)
        elseif length(ids) == 4
            # Preserve QUAD
            push!(cell_node_ids, ids)
            push!(cell_types, 2) # 2 = QUAD
            current_cell_idx += 1
            push!(created_cell_indices, current_cell_idx)
        else
            error("[GridapBridge] Unsupported element type with $(length(ids)) nodes. Only TRI and QUAD supported.")
        end

        # Record Lineage (Accumulate, do not overwrite!)
        if !haskey(leaf_to_cell_ids, el.leaf_id)
            leaf_to_cell_ids[el.leaf_id] = Int[]
        end
        append!(leaf_to_cell_ids[el.leaf_id], created_cell_indices)
    end

    # 3. Construct Gridap Data Structures
    # -----------------------------------
    
    n_cells = length(cell_node_ids)
    n_nodes = length(node_coords)
    
    data = Int32[]
    ptrs = Int32[1]
    
    for ids in cell_node_ids
        append!(data, Int32.(ids))
        push!(ptrs, ptrs[end] + Int32(length(ids)))
    end
    cell_to_nodes_table = Table(data, ptrs)
    
    # Helper for Reffes
    # We define the available types dynamically
    
    reffe_tri = LagrangianRefFE(Float64, TRI, 1)
    reffe_quad = LagrangianRefFE(Float64, QUAD, 1)
    
    has_tri = any(length(ids) == 3 for ids in cell_node_ids)
    has_quad = any(length(ids) == 4 for ids in cell_node_ids)
    
    unique_reffes = Vector{LagrangianRefFE{2}}()
    
    # Map from local type (1=Tri, 2=Quad in loop above) to Index in unique_reffes
    # We used: 1 -> TRI, 2 -> QUAD
    type_map = Dict{Int, Int}()
    
    if has_tri
        push!(unique_reffes, reffe_tri)
        type_map[1] = length(unique_reffes)
    end
    
    if has_quad
        push!(unique_reffes, reffe_quad)
        type_map[2] = length(unique_reffes)
    end
    
    # Remap indices and force Int8
    cell_type_indices = Int8[type_map[t] for t in cell_types]
    
    # Create Unstructured Grid
    # Safety Checks
    if isempty(cell_types)
        error("[GridapBridge] No elements provided!")
    end
    @assert length(cell_types) == length(cell_node_ids)
    
    # Check for degenerate elements
    for (i, ids) in enumerate(cell_node_ids)
        if any(x -> x <= 0, ids)
            error("[GridapBridge] Cell $i has invalid node IDs: $ids")
        end
        if length(unique(ids)) != length(ids)
             error("[GridapBridge] Degenerate Cell $i (duplicate nodes): $ids. Coords: $([node_coords[id] for id in ids])")
        end
    end
    
    println("  [GridapBridge] Element Types Summary: TRI=$(count(==(1), cell_types)), QUAD=$(count(==(2), cell_types))")


    grid = UnstructuredGrid(
        node_coords,
        cell_to_nodes_table,
        unique_reffes,
        cell_type_indices,
        Gridap.Geometry.NonOriented()
    )
    
    # Create Model
    model = UnstructuredDiscreteModel(grid)
    
    println("  [GridapBridge] Node Merging: Queries=$n_queries, Merged=$n_merged, Ratio=$(n_merged/n_queries)")
    println("  [GridapBridge] Success! Model has $(num_cells(model)) cells (Mixed TRI/QUAD) and $(num_nodes(model)) nodes.")
    return model, leaf_to_cell_ids
end

end # module
