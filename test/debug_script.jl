using Gridap
using Gridap.Geometry

# Minimal QuadNode stub
mutable struct QuadNode
    id::Int
    center::Vector{Float64}
    size::Float64
    level::Int
    is_active::Bool
    children::Vector{QuadNode}
    parent::Union{QuadNode, Nothing}
    QuadNode(id, center, size, level) = new(id, center, size, level, true, QuadNode[], nothing)
end

function run_debug()
    domain = (0, 1, 0, 1)
    partition = (50, 50)
    model = CartesianDiscreteModel(domain, partition)
    
    desc = get_cartesian_descriptor(model)
    origin = desc.origin
    sizes = desc.sizes
    partition = desc.partition
    
    nx, ny = partition[1], partition[2]
    Lx, Ly = sizes[1], sizes[2]
    hx = Lx / nx
    hy = Ly / ny
    cell_size = (hx + hy) / 2.0
    
    println("Lx=$Lx, nx=$nx, hx=$hx, cell_size=$cell_size")
    
    # Simulate Loop
    current_depth = 0
    current_nx = nx
    current_ny = ny
    current_size = cell_size
    
    while current_nx > 1 || current_ny > 1
        parent_nx = ceil(Int, current_nx / 2)
        parent_ny = ceil(Int, current_ny / 2)
        
        println("Depth=$current_depth. nx=$current_nx. current_size=$current_size. p_size expected=$(current_size*2.0)")
        
        current_depth += 1
        current_nx = parent_nx
        current_ny = parent_ny
        current_size *= 2.0
    end
    println("Final Size: $current_size")
end

run_debug()
