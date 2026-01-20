using CoarseAgFEM
using Gridap
using Gridap.Geometry
using GridapEmbedded
using Test

function test_crescent()
    println("=== Test: Crescent Shape (500x500) ===")

    # 1. Create Input Cartesian Model
    # --------------------------------
    # High resolution: 500x500
    domain = (0, 1, 0, 1)
    partition = (500, 500) 
    model = CartesianDiscreteModel(domain, partition)
    
    println("  Input Model: $(num_cells(model)) cells.")
    
    # 2. Convert to QuadMesh
    # ----------------------
    qmesh = cartesian_to_quadtree(model)
    
    # Check leaves count
    leaves = [n for n in qmesh.all_nodes if n.is_active && isempty(n.children)]
    println("  QuadMesh Leaves: $(length(leaves))")
    @test length(leaves) == 250000

    # 3. Apply Coarsening (Greedy Packing)
    # ------------------------------------
    # Crescent Level Set
    # Intersection of Circle A (In) and Circle B (Out)
    # Phi = max(Phi_A, -Phi_B)
    
    function level_set_func(x)
        cA = (0.5, 0.5); rA = 0.4
        cB = (0.6, 0.5); rB = 0.3
        
        dA = sqrt((x[1]-cA[1])^2 + (x[2]-cA[2])^2)
        dB = sqrt((x[1]-cB[1])^2 + (x[2]-cB[2])^2)
        
        phi_A = dA - rA
        phi_B = dB - rB
        
        return max(phi_A, -phi_B)
    end
    
    # Step 3a: Classify
    # Use 4-cell buffer safety margin (h=0.002)
    h_min = 0.002
    buffer_margin = 4.0 * h_min
    classify_leaves!(qmesh, level_set_func; buffer_width=buffer_margin)
    
    # Step 3b: Greedy Homogeneous Coarsening
    # Limit max cell size to 8 * h_min
    bottom_up_coarsening!(qmesh; max_coarsening_factor=8)
    
    new_leaves = [n for n in qmesh.all_nodes if n.is_active && isempty(n.children)]
    println("  Coarsened Leaves: $(length(new_leaves))")
    
    # Debug: Count by status
    _INTERIOR = CoarseAgFEM.INTERIOR
    _EXTERIOR = CoarseAgFEM.EXTERIOR
    _CUT      = CoarseAgFEM.CUT
    _BUFFER   = CoarseAgFEM.BUFFER
    
    c_int = count(n -> n.status == _INTERIOR, new_leaves)
    c_ext = count(n -> n.status == _EXTERIOR, new_leaves)
    c_cut = count(n -> n.status == _CUT, new_leaves)
    c_buf = count(n -> n.status == _BUFFER, new_leaves)
    println("  Leaves Stats: Int=$c_int, Ext=$c_ext, Cut=$c_cut, Buf=$c_buf")
    
    # 4. Balance
    # ----------
    balance!(qmesh)
    
    # 5. Connect Paving (Skipped)
    # ---------------------------
    println("  Phase 3: Paving SKIPPED (User Request). Generating raw Quad elements...")
    
    active_leaves = [n for n in qmesh.all_nodes if n.is_active && isempty(n.children)]
    elements = QuadElement[]
    for leaf in active_leaves
        b = get_bounds(leaf)
        # Counter-Clockwise Order: SW, SE, NE, NW
        p1 = [b[1], b[3]]
        p2 = [b[2], b[3]]
        p3 = [b[2], b[4]]
        p4 = [b[1], b[4]]
        push!(elements, QuadElement(leaf.id, [p1, p2, p3, p4], "orange"))
    end
    
    println("  Generated $(length(elements)) elements (All Quads).")
    
    # 6. Convert to Gridap Model (Validity Check)
    # -------------------------------------------
    out_model, lineage = quadtree_to_discrete_model(elements)
    println("  Output Model: $(num_cells(out_model)) cells.")
    
    # Write VTK (Internal)
    mkpath("output")
    # write_vtk("output/crescent_test_mesh.vtu", elements, qmesh)
    writevtk(out_model, "output/crescent_test_mesh")
    println("  Wrote output/crescent_test_mesh.vtu")
    
    @test num_cells(out_model) > 0
end

test_crescent()
