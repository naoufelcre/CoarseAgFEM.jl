using CoarseAgFEM
using Gridap
using Gridap.Geometry
using GridapEmbedded
using Gridap.ReferenceFEs
using Test

function compare_orphans()
    println("=== Comparison: Uniform vs Quadtree Aggregation Robustness (Circle) ===")
    
    cases = [
        ("Circle", x -> sqrt((x[1]-0.5)^2 + (x[2]-0.5)^2) - 0.4),
        ("Crescent", x -> begin
            cA = (0.5, 0.5); rA = 0.4
            cB = (0.6, 0.5); rB = 0.3
            dA = sqrt((x[1]-cA[1])^2 + (x[2]-cA[2])^2)
            dB = sqrt((x[1]-cB[1])^2 + (x[2]-cB[2])^2)
            max(dA - rA, -(dB - rB))
        end)
    ]

    for (name, ls_func) in cases
        println("\n=== Comparison: $name ===")
        geo = AnalyticalGeometry(x -> ls_func(x))
        
        # ---------------------------------------------------------
        # Case A: Uniform 500x500
        # ---------------------------------------------------------
        println("  [Case A] Uniform 500x500")
        domain = (0, 1, 0, 1)
        partition = (500, 500)
        model_uni = CartesianDiscreteModel(domain, partition)
        
        cutgeo_uni = cut(model_uni, geo)
        
        strategy = RobustAggregation(0.5, 1000)
        aggs_uni = aggregate(strategy, cutgeo_uni, geo, IN)
        
        # Orphans
        cell_to_inoutcut_uni = GridapEmbedded.AgFEM.compute_bgcell_to_inoutcut(cutgeo_uni, geo)
        _CUT_val = GridapEmbedded.Interfaces.CUT
        orphans_uni = count(i -> (aggs_uni[i] == 0 && cell_to_inoutcut_uni[i] == _CUT_val), 1:length(aggs_uni))
        println("    -> Orphans: $orphans_uni")
        
        # Output
        # Improved Visualization: Separate Orphans from Exterior
        # aggs[i] == 0 for Exterior cells too! We must distinguish.
        
        is_cut_orphan_uni = zeros(Float64, length(aggs_uni))
        for i in 1:length(aggs_uni)
            if aggs_uni[i] == 0 && cell_to_inoutcut_uni[i] == _CUT_val
                is_cut_orphan_uni[i] = 1.0
            end
        end
        
        fields_uni = Dict(
            "aggs" => aggs_uni, 
            "inoutcut" => cell_to_inoutcut_uni,
            "is_cut_orphan" => is_cut_orphan_uni
        )
        
        mkpath("output")
        writevtk(Triangulation(model_uni), "output/$(name)_uniform_aggregation"; celldata=fields_uni)
        
        try
            trian_phys = Triangulation(cutgeo_uni)
            writevtk(trian_phys, "output/$(name)_uniform_physical")
        catch e
             println("    [Warn] Failed to write physical domain: $e")
        end

        # ---------------------------------------------------------
        # Case B: Quadtree 8x
        # ---------------------------------------------------------
        println("  [Case B] Quadtree 8x")
        qmesh = initialize_builder(model_uni)
        
        h_min = 0.002
        buffer_margin = 4.0 * h_min
        classify_leaves!(qmesh, ls_func; buffer_width=buffer_margin)
        bottom_up_coarsening!(qmesh; max_coarsening_factor=8)
        balance!(qmesh)
        
        active_leaves = [n for n in qmesh.all_nodes if n.is_active && isempty(n.children)]
        elements = QuadElement[]
        for leaf in active_leaves
            b = CoarseAgFEM.QuadtreeMeshing.get_bounds(leaf)
            p1 = [b[1], b[3]]; p2 = [b[2], b[3]]; p3 = [b[2], b[4]]; p4 = [b[1], b[4]]
            # Gridap QUAD Node Ordering: 1=BL, 2=BR, 3=TL, 4=TR !!
            # We defined p3=TR, p4=TL. So we must push [p1, p2, p4, p3].
            push!(elements, QuadElement(leaf.id, [p1, p2, p4, p3], "orange"))
        end
        
        model_qt, _ = quadtree_to_discrete_model(elements)
        
        cutgeo_qt = cut(model_qt, geo)
        aggs_qt = aggregate(strategy, cutgeo_qt, geo, IN)
        
        cell_to_inoutcut_qt = GridapEmbedded.AgFEM.compute_bgcell_to_inoutcut(cutgeo_qt, geo)
        orphans_qt = count(i -> (aggs_qt[i] == 0 && cell_to_inoutcut_qt[i] == _CUT_val), 1:length(aggs_qt))
        println("    -> Orphans: $orphans_qt")
        
        # --- CRITICAL CHECK: Are any COARSE cells cut? ---
        # If a coarsened cell (Level > 0) is CUT, our classification failed.
        # We can check cell volumes.
        trian_qt = Triangulation(model_qt)
        measures = get_cell_measure(trian_qt)
        # Expected fine area = (0.002)^2 = 4.0e-6
        # Allow small epsilon.
        fine_area = 4.0e-6
        tol = 1.0e-8
        
        coarse_cuts = 0
        for i in 1:length(cell_to_inoutcut_qt)
            if cell_to_inoutcut_qt[i] == _CUT_val
                if measures[i] > (fine_area + tol)
                    coarse_cuts += 1
                end
            end
        end
        
        if coarse_cuts > 0
            println("    [CRITICAL FAIL] Found $coarse_cuts COARSE cells that are CUT! Buffer Zone Violation.")
        else
            println("    [PASS] All CUT cells are at finest level (Buffer Zone respected).")
        end
        # -------------------------------------------------
        
        is_cut_orphan_qt = zeros(Float64, length(aggs_qt))
        for i in 1:length(aggs_qt)
            if aggs_qt[i] == 0 && cell_to_inoutcut_qt[i] == _CUT_val
                is_cut_orphan_qt[i] = 1.0
            end
        end
        
        fields_qt = Dict(
            "aggs" => aggs_qt, 
            "inoutcut" => cell_to_inoutcut_qt, 
            "is_cut_orphan" => is_cut_orphan_qt
        )
        
        mkpath("output")
        writevtk(Triangulation(model_qt), "output/$(name)_quadtree_aggregation"; celldata=fields_qt)
        
        try
            trian_phys_qt = Triangulation(cutgeo_qt)
            writevtk(trian_phys_qt, "output/$(name)_quadtree_physical")
        catch e
            println("    [Warn] Failed to write physical domain: $e")
        end
    end
end

compare_orphans()

compare_orphans()
