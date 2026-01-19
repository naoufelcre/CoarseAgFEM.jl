# Crescent Economy Example
# This script demonstrates the "Cell Economy" of CoarseAgFEM.
# It compares a high-resolution Uniform Mesh vs. an Adaptive Quadtree Mesh.

using CoarseAgFEM
using Gridap
using Gridap.Geometry
using Printf

function crescent_economy_example()
    println("==========================================")
    println("   CoarseAgFEM: Crescent Economy Example  ")
    println("==========================================")

    mkpath("output")

    # 1. Geometry Definition (Crescent)
    # ---------------------------------
    # Intersection of Circle A (Inside) and Circle B (Outside)
    # Phi = max(Phi_A, -Phi_B)
    cA = (0.5, 0.5); rA = 0.4
    cB = (0.7, 0.5); rB = 0.3
    
    function level_set_func(x)
        dA = sqrt((x[1]-cA[1])^2 + (x[2]-cA[2])^2)
        dB = sqrt((x[1]-cB[1])^2 + (x[2]-cB[2])^2)
        phi_A = dA - rA
        phi_B = dB - rB
        return max(phi_A, -phi_B)
    end

    # 2. Parameters
    # -------------
    # We want a fine resolution at the interface (h_min)
    # Target Resolution: 1000x1000 grid equivalent
    n_fine = 1000
    domain = (0, 1, 0, 1)
    
    println("\n--- 1. Uniform Cartesian Mesh (Reference) ---")
    model_uniform = CartesianDiscreteModel(domain, (n_fine, n_fine))
    n_cells_uniform = num_cells(model_uniform)
    println("  Grid Resolution: $n_fine x $n_fine")
    println("  Total Cells: $n_cells_uniform")
    
    # 3. Quadtree Mesh (Adaptive)
    # ---------------------------
    println("\n--- 2. Adaptive Quadtree Mesh ---")
    
    # Start with the same fine grid as leaves
    println("  Initializing Quadtree from Uniform Grid...")
    qmesh = cartesian_to_quadtree(model_uniform)
    
    # A. Classification
    # Use a small buffer to ensure interface capture
    h_min = 1.0 / n_fine
    buffer_width = 4.0 * h_min
    println("  Classifying leaves (Buffer: $(buffer_width))...")
    classify_leaves!(qmesh, level_set_func; buffer_width=buffer_width)
    
    # B. Coarsening
    # Allow coarsening away from interface up to factor 32 (Level 5)
    max_factor = 32
    println("  Coarsening (Max Factor: $max_factor)...")
    bottom_up_coarsening!(qmesh; max_coarsening_factor=max_factor)
    
    # C. Balancing
    println("  Balancing mesh (2:1 rule)...")
    balance!(qmesh)
    
    # D. Paving
    println("  Paving elements...")
    elements = pave_mesh(qmesh)
    
    # E. Convert to Gridap
    model_quadtree, _ = quadtree_to_discrete_model(elements)
    n_cells_quadtree = num_cells(model_quadtree)
    
    println("  Total Cells: $n_cells_quadtree")

    # 4. Economy Analysis
    # -------------------
    ratio = n_cells_uniform / n_cells_quadtree
    reduction = 100.0 * (1.0 - n_cells_quadtree / n_cells_uniform)
    
    println("\n==========================================")
    println("             RESULTS SUMMARY              ")
    println("==========================================")
    @printf "  Uniform Cells:  %d\n" n_cells_uniform
    @printf "  Quadtree Cells: %d\n" n_cells_quadtree
    @printf "  Reduction:      %.2f%%\n" reduction
    @printf "  Speedup Factor: %.2fx (Theoretical)\n" ratio
    println("==========================================\n")

    # 5. Visualization
    # ----------------
    println("  Writing visualizations to output/...")
    try
        writevtk(model_uniform, "output/crescent_uniform")
        println("  -> Written output/crescent_uniform.vtu")
    catch e
        println("  [Warning] Failed to write Uniform VTK: $e")
    end

    try
        writevtk(model_quadtree, "output/crescent_quadtree")
        println("  -> Written output/crescent_quadtree.vtu")
    catch e
        println("  [Warning] Failed to write Quadtree VTK (Gridap): $e")
        println("  -> Attempting fallback internal writer...")
        try
            CoarseAgFEM.write_vtk("output/crescent_quadtree_raw.vtu", elements, qmesh)
            println("  -> Written output/crescent_quadtree_raw.vtu (Fallback)")
        catch e2
            println("  [Error] Fallback writer also failed: $e2")
        end
    end
end

crescent_economy_example()
