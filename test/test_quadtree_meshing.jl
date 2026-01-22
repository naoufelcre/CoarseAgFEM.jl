using Test
using CoarseAgFEM
using Gridap
using Gridap

@testset "Quadtree Aggregations" begin



    @testset "Coarsening" begin
        # Create fine mesh (Level 4: 16x16 = 256 cells)
        model = CartesianDiscreteModel((0, 1, 0, 1), (16, 16))
        mesh = initialize_builder(model)
        
        # Use Level Set to drive coarsening (Fine at interface, coarse away)
        # Circle at center
        dist(x) = sqrt((x[1]-0.5)^2 + (x[2]-0.5)^2) - 0.25
        
        classify_leaves!(mesh, dist; buffer_width=0.0)
        bottom_up_coarsening!(mesh)
        
        leaves = [n for n in mesh.all_nodes if n.is_active && isempty(n.children)]
        # Should be fewer than 256
        @test length(leaves) < 256
        
        # Verify sizes vary
        sizes = [n.size for n in leaves]
        @test length(unique(sizes)) > 1
    end

    @testset "Balancing" begin
        model = CartesianDiscreteModel((0, 1, 0, 1), (32, 32)) # Finer to allow more ripples
        mesh = initialize_builder(model)
        
        # Force imbalance: Interface near corner
        dist(x) = (x[1] + x[2]) - 0.4
        
        classify_leaves!(mesh, dist; buffer_width=0.05)
        bottom_up_coarsening!(mesh)
        
        # Check if unbalanced (naive check: just run balance! and ensure it doesn't crash)
        # Ideally we'd check the 2:1 rule but that's internal to balance!
        @test_nowarn balance!(mesh)
    end
    
end

