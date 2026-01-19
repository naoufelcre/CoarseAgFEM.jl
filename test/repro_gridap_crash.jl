using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.Arrays

function test_mixed_mesh_vtk()
    println("Testing Mixed Mesh VTK generation...")

    # Define nodes: square (0,0) to (1,1)
    # 1 -- 2
    # |  / |
    # 3 -- 4
    #
    # Triangle: 1, 3, 2
    # Quad: (using 2,3,4? No, let's make a separate Quad)
    #
    # Let's do:
    # Tri: (0,0), (1,0), (0,1) -> Nodes 1, 2, 3
    # Quad: (1,0), (2,0), (2,1), (1,1) -> Nodes 2, 4, 5, 6
    
    node_coords = [
        Point(0.0, 0.0), # 1
        Point(1.0, 0.0), # 2
        Point(0.0, 1.0), # 3
        Point(2.0, 0.0), # 4
        Point(2.0, 1.0), # 5
        Point(1.0, 1.0)  # 6
    ]
    
    # Cells
    # 1: TRI [1, 2, 3]
    # 2: QUAD [2, 4, 5, 6]
    
    data = [1, 2, 3, 2, 4, 5, 6]
    ptrs = [1, 4, 8]
    cell_to_nodes = Table(data, ptrs)
    
    reffe_tri = LagrangianRefFE(Float64, TRI, 1)
    reffe_quad = LagrangianRefFE(Float64, QUAD, 1)
    
    # Note: explicit typing to match GridapIntegration
    unique_reffes = Vector{LagrangianRefFE{2}}([reffe_tri, reffe_quad])
    
    # Cell 1 -> Type 1 (Tri)
    # Cell 2 -> Type 2 (Quad)
    cell_types = Int8[1, 2]
    
    grid = UnstructuredGrid(
        node_coords,
        cell_to_nodes,
        unique_reffes,
        cell_types,
        Gridap.Geometry.NonOriented()
    )
    
    model = UnstructuredDiscreteModel(grid)
    
    mkpath("output")
    println("Writing VTK...")
    writevtk(model, "output/repro_mixed_mesh")
    println("Success!")
end

test_mixed_mesh_vtk()
