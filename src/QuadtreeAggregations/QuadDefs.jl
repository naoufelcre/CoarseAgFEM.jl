module QuadDefs

export QuadNode, QuadMesh, QuadElement
export is_leaf, get_bounds

@enum NodeStatus INTERIOR EXTERIOR CUT BUFFER UNDEFINED

mutable struct QuadNode
    id::Int
    parent::Union{QuadNode, Nothing}
    children::Vector{QuadNode} 
    
    center::Vector{Float64}
    size::Float64
    level::Int
    
    # helper for fast merging
    is_active::Bool 
    status::NodeStatus
    
    QuadNode(id, center, size, level) = new(id, nothing, QuadNode[], center, size, level, true, UNDEFINED)
end

export NodeStatus, INTERIOR, EXTERIOR, CUT, BUFFER, UNDEFINED

# Renamed to avoid specific conflict, though Element is generic.
struct QuadElement
    leaf_id::Int  # Lineage: Which QuadNode created this?
    nodes::Vector{Vector{Float64}}
    color::String 
end

mutable struct QuadMesh
    root::QuadNode
    all_nodes::Vector{QuadNode} # Flat list of ALL nodes created (including deleted ones)
end

is_leaf(n::QuadNode) = isempty(n.children)

function get_bounds(n::QuadNode)
    h = n.size / 2.0
    return (n.center[1]-h, n.center[1]+h, n.center[2]-h, n.center[2]+h)
end

end # module
