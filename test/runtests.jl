using Test
using CoarseAgFEM

@testset "CoarseAgFEM" begin
    include("test_robust_agfem.jl")
    include("test_quadtree_meshing.jl")
    include("test_refactor.jl")
    include("test_aggregation.jl")
    include("test_crescent.jl")
end
