using QuadtreeAgFEM
using Gridap
using Gridap.Geometry
using GridapEmbedded
using Test

function test_refactor()
    println("=== Test: Uniform-First Quadtree Refactor ===")

    # 1. Create Input Cartesian Model
    # --------------------------------
    # Non-POT grid: 500x500
    domain = (0, 1, 0, 1)
    partition = (500, 500) 
    model = CartesianDiscreteModel(domain, partition)
    
    println("  Input Model: $(num_cells(model)) cells.")

    # 2. Convert to QuadMesh
    # ----------------------
    qmesh = cartesian_to_quadtree(model)
    
    # Check leaves count
    leaves = [n for n in qmesh.all_nodes if n.is_active && isempty(n.children)]
    println("  QuadMesh Leaves: $(length(leaves)) (Expected 250000)")
    @test length(leaves) == 250000

    # 3. Apply Coarsening (Greedy Packing)
    # ------------------------------------
    # Circle Level Set: Center (0.5, 0.5), Radius 0.2
    # Negative inside, Positive outside.
    
    function level_set_func(x)
        # distance - radius
        return sqrt((x[1]-0.5)^2 + (x[2]-0.5)^2) - 0.2
    end
    
    # Step 3a: Classify
    # Use 4-cell buffer safety margin
    h_min = 1.0 / 500.0 # 0.002
    buffer_margin = 4.0 * h_min
    classify_leaves!(qmesh, level_set_func; buffer_width=buffer_margin)
    
    # Step 3b: Greedy Homogeneous Coarsening
    bottom_up_coarsening!(qmesh; max_coarsening_factor=8)
    
    new_leaves = [n for n in qmesh.all_nodes if n.is_active && isempty(n.children)]
    println("  Coarsened Leaves: $(length(new_leaves))")
    
    # Debug: Count by status
    # Debug: Count by status
    # Access via internal path to ensure availability
    _INTERIOR = QuadtreeAgFEM.QuadtreeAggregations.QuadDefs.INTERIOR
    _EXTERIOR = QuadtreeAgFEM.QuadtreeAggregations.QuadDefs.EXTERIOR
    _CUT      = QuadtreeAgFEM.QuadtreeAggregations.QuadDefs.CUT
    _BUFFER   = QuadtreeAgFEM.QuadtreeAggregations.QuadDefs.BUFFER
    
    c_int = count(n -> n.status == _INTERIOR, new_leaves)
    c_ext = count(n -> n.status == _EXTERIOR, new_leaves)
    c_cut = count(n -> n.status == _CUT, new_leaves)
    c_buf = count(n -> n.status == _BUFFER, new_leaves)
    println("  Leaves Stats: Int=$c_int, Ext=$c_ext, Cut=$c_cut, Buf=$c_buf")
    
    # 4. Balance (Optional but good for paving)
    # ----------------------------------------
    balance!(qmesh)
    
    # 5. Connect Paving (Skipped)
    # ---------------------------
    println("  Phase 3: Paving SKIPPED (User Request). Generating raw Quad elements...")
    
    active_leaves = [n for n in qmesh.all_nodes if n.is_active && isempty(n.children)]
    elements = QuadElement[]
    for leaf in active_leaves
        b = QuadtreeAgFEM.QuadtreeAggregations.get_bounds(leaf) # Deep Qualified or use exported
        # Counter-Clockwise Order: SW, SE, NE, NW
        p1 = [b[1], b[3]]
        p2 = [b[2], b[3]]
        p3 = [b[2], b[4]]
        p4 = [b[1], b[4]]
        push!(elements, QuadElement(leaf.id, [p1, p2, p3, p4], "orange"))
    end
    
    println("  Generated $(length(elements)) elements (All Quads).")
    
    # 6. Convert to Gridap Model
    # --------------------------
    out_model, lineage = quadtree_to_discrete_model(elements)
    
    # Verify outputs
    println("  Output Model: $(num_cells(out_model)) cells.")
    
    # Check if we have both QUAD (dim 2, type ?) and TRI 
    # Gridap doesn't easily expose cell types summary without iteration
    # But we can write VTK
    
    writevtk(out_model, "output/gridap_refactor_test_mesh")
    # Gridap vtk crashes on mixed mesh currently.
    
    # Use internal VTK writer which is robust for QuadMesh
    # Signature: write_vtk(filename, elements, mesh)
    # QuadtreeAgFEM.write_vtk("output/refactor_test_mesh.vtu", elements, qmesh)
    println("  Wrote output/gridap_refactor_test_mesh.vtu")
    
    # Basic Assertion: We didn't crash.
    @test num_cells(out_model) > 0
end

test_refactor()
