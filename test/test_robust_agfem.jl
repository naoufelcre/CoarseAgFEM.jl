using Test
using CoarseAgFEM
using Gridap
using GridapEmbedded
using GridapEmbedded.LevelSetCutters

@testset "Robust AgFEM" begin
    # Create simple Cartesian model
    n = 20
    model = CartesianDiscreteModel((0,1,0,1), (n,n))
    
    # Define a cut that creates some small elements (circle)
    R = 0.33452 # Non-rational radius to avoid perfect symmetry
    circle_sdf(x) = sqrt((x[1]-0.5)^2 + (x[2]-0.5)^2) - R
    
    # Cut
    geo = DiscreteGeometry(circle_sdf, model)
    cut_geo = cut(model, geo)
    
    # 1. Standard Aggregation (Baseline)
    strategy_std = AggregateCutCellsByThreshold(0.5)
    aggs_std = aggregate(strategy_std, cut_geo)
    @test length(aggs_std) == num_cells(model)
    
    # 2. Robust Aggregation
    strategy_robust = RobustAggregation(0.5, 50) # Explicit higher limit
    aggs_robust = aggregate(strategy_robust, cut_geo)
    
    @test length(aggs_robust) == num_cells(model)
    @test aggs_robust == aggs_std # Should be identical for simple case
    
    println("  Robust AgFEM Test Passed (Standard Symmetry)")
end
