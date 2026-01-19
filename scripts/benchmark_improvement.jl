using QuadtreeAgFEM
using Gridap
using Printf
using LinearAlgebra

# ----------------------------------------------------------------------
# Helper Functions
# ----------------------------------------------------------------------

function circle_sdf(x)
    return norm(x .- [0.5,0.5]) - 0.35
end

function crescent_sdf(x)
    c1 = norm(x .- [0.5, 0.5]) - 0.4
    c2 = norm(x .- [0.6, 0.5]) - 0.3
    return max(c1, -c2)
end

function run_benchmark(N::Int, geometry_name::String, dist_func::Function; 
                       max_factor::Int=16, 
                       buffer_cells::Int=3)
    
    println("\n======================================================================")
    println("BENCHMARK: $geometry_name ($N x $N)")
    println("Params: Max Coarsening Factors=$max_factor, Buffer Cells=$buffer_cells")
    println("======================================================================")
    
    h_fine = 1.0 / N
    buf_width = buffer_cells * h_fine
    
    # 0. Fine Grid Generation (Implicit)
    t0 = time()
    cmodel = CartesianDiscreteModel((0,1,0,1), (N,N))
    mesh = cartesian_to_quadtree(cmodel)
    t_gen = time() - t0
    println("  Generation (Fine):    $(Printf.@sprintf("%.4f", t_gen)) s")
    
    # 1. Classification
    t1 = time()
    classify_leaves!(mesh, dist_func, buffer_width=buf_width)
    
    # 2. Coarsening
    bottom_up_coarsening!(mesh, max_coarsening_factor=max_factor)
    t_coarse = time() - t1
    println("  Coarsening:           $(Printf.@sprintf("%.4f", t_coarse)) s")
    
    # 3. Balancing
    t2 = time()
    balance!(mesh)
    t_bal = time() - t2
    println("  Balancing:            $(Printf.@sprintf("%.4f", t_bal)) s")
    
    # 4. Paving
    t3 = time()
    elements = pave_mesh(mesh)
    t_pave = time() - t3
    println("  Paving (Hybrid):      $(Printf.@sprintf("%.4f", t_pave)) s")
    
    # 5. Gridap Conversion
    t4 = time()
    model, lineage = quadtree_to_discrete_model(elements)
    t_conv = time() - t4
    println("  Gridap Conversion:    $(Printf.@sprintf("%.4f", t_conv)) s")
    
    # Stats
    n_uniform = N*N
    n_quads = num_cells(model)
    reduction = n_uniform / n_quads
    total_time = t_gen + t_coarse + t_bal + t_pave + t_conv
    
    println("\n  RESULTS:")
    println("  - Uniform Cells:      $n_uniform")
    println("  - Quadtree Cells:     $n_quads")
    println("  - Reduction Factor:   $(Printf.@sprintf("%.2f", reduction)) x ($(Printf.@sprintf("%.2f", 100*n_quads/n_uniform))% of original)")
    println("  - Total Time:         $(Printf.@sprintf("%.4f", total_time)) s")
    
    # Verify Topology
    grid = get_grid(model)
    cell_node_ids = Gridap.Geometry.get_cell_node_ids(grid)
    n_tris = count(ids -> length(ids) == 3, cell_node_ids)
    
    println("  - Analysis:           $(n_quads - n_tris) Quads, $n_tris Triangles")
    
    return n_quads
end

# ----------------------------------------------------------------------
# Execution (High Resolution)
# ----------------------------------------------------------------------

# Warmup
# run_benchmark(64, "Warmup", circle_sdf)

# High Res
# run_benchmark(1024, "Circle", circle_sdf, max_factor=16, buffer_cells=3)
run_benchmark(1024, "Crescent", crescent_sdf, max_factor=16, buffer_cells=3)
run_benchmark(2048, "Crescent", crescent_sdf, max_factor=16, buffer_cells=3)

