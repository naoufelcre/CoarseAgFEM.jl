# Crescent Economy Example
# This script demonstrates the "Cell Economy" of CoarseAgFEM.
# It compares a high-resolution Uniform Mesh vs. an Adaptive Quadtree Mesh.

using CoarseAgFEM
using Gridap
using Gridap.Geometry
using Printf

function crescent_economy_example(n_fine::Int)
    println("\n==========================================")
    println("   CoarseAgFEM: Crescent @ $(n_fine)x$(n_fine)")
    println("==========================================")

    mkpath("output")

    # 1. Geometry Definition (Crescent)
    # ---------------------------------
    cA = (0.5, 0.5); rA = 0.4
    cB = (0.7, 0.5); rB = 0.3
    
    function level_set_func(x)
        dA = sqrt((x[1]-cA[1])^2 + (x[2]-cA[2])^2)
        dB = sqrt((x[1]-cB[1])^2 + (x[2]-cB[2])^2)
        phi_A = dA - rA
        phi_B = dB - rB
        return max(phi_A, -phi_B)
    end

    domain = (0, 1, 0, 1)
    
    println("\n--- 1. Uniform Cartesian Mesh (Reference) ---")
    model_uniform = CartesianDiscreteModel(domain, (n_fine, n_fine))
    n_cells_uniform = num_cells(model_uniform)
    println("  Grid Resolution: $n_fine x $n_fine")
    println("  Total Cells: $n_cells_uniform")
    
    # 2. Quadtree Mesh (Adaptive)
    # ---------------------------
    println("\n--- 2. Adaptive Quadtree Mesh ---")
    
    h_min = 1.0 / n_fine
    buffer_width = 4.0 * h_min
    max_factor = 32
    
    builder = initialize_builder(model_uniform)
    classify_leaves!(builder, level_set_func; buffer_width=buffer_width)
    bottom_up_coarsening!(builder; max_coarsening_factor=max_factor)
    balance!(builder)
    
    elements = pave_mesh(builder)
    model_quadtree, _ = quadtree_to_discrete_model(elements)
                                            
    n_cells_quadtree = num_cells(model_quadtree)
    println("  Total Cells: $n_cells_quadtree")

    # 3. Economy Analysis
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

    # 4. Visualization
    # ----------------
    println("  Writing SVG to output/crescent_$(n_fine).svg ...")
    write_svg(model_quadtree, "output/crescent_$(n_fine)")
    println("  -> Done!")
    
    return (n_fine, n_cells_uniform, n_cells_quadtree, reduction, ratio)
end

# Run for multiple scales
results = []
for scale in [250, 500, 1000]
    push!(results, crescent_economy_example(scale))
end

# Summary Table
println("\n==========================================")
println("         MULTI-SCALE SUMMARY              ")
println("==========================================")
println("  Scale   | Uniform   | Quadtree | Reduction | Speedup")
println("  --------|-----------|----------|-----------|--------")
for (n, uni, qt, red, spd) in results
    @printf "  %4dx%4d | %9d | %8d | %6.2f%%   | %5.2fx\n" n n uni qt red spd
end
println("==========================================")

