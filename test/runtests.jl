using Test
using CoarseAgFEM

@testset "CoarseAgFEM" begin
    include("test_robust_agfem.jl")
    include("test_quadtree_meshing.jl")
end
