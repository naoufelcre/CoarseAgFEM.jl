using CoarseAgFEM
using Gridap

# 1. Create a simple model
domain = (0, 1, 0, 1)
partition = (4, 4)
model = CartesianDiscreteModel(domain, partition)

# 2. Build coarse model
builder = initialize_builder(model)
elements = pave_mesh(builder)

# 3. Test write_svg (if exported or via QuadtreeMeshing)
# Since it's not exported from CoarseAgFEM, we test if we can access it
try
    CoarseAgFEM.write_svg("test_output.svg", elements)
    println("SUCCESS: write_svg worked via CoarseAgFEM (wait, it shouldn't be exported?)")
catch e
    println("EXPECTED: write_svg not directly in CoarseAgFEM: $e")
    # Try via QuadtreeMeshing
    CoarseAgFEM.QuadtreeMeshing.write_svg("test_output.svg", elements)
    println("SUCCESS: write_svg worked via CoarseAgFEM.QuadtreeMeshing")
end
