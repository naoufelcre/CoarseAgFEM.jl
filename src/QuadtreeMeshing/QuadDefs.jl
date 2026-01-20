module QuadDefs

export QuadNode, CoarseMeshBuilder, QuadElement
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

"""
    QuadNode (Builder Node)

Represents a node in the construction tree/forest.
Used only by `CoarseMeshBuilder`.
"""

export NodeStatus, INTERIOR, EXTERIOR, CUT, BUFFER, UNDEFINED

struct QuadElement
    leaf_id::Int  # Lineage: Which QuadNode created this?
    nodes::Vector{Vector{Float64}}
    color::String
end

"""
    mutable struct CoarseMeshBuilder

Transient builder object for constructing a coarse mesh hierarchy.
Stores the forest of QuadNodes and manages the coarsening/balancing process.
NOT INTENDED TO BE USED AS A PERMANENT MESH STRUCTURE.
"""
mutable struct CoarseMeshBuilder
    roots::Vector{QuadNode}     # Support for "Forest" of roots (non-POT domains)
    all_nodes::Vector{QuadNode} # Flat list of ALL nodes created (including deleted ones)
end

is_leaf(n::QuadNode) = isempty(n.children)

function get_bounds(n::QuadNode)
    h = n.size / 2.0
    return (n.center[1]-h, n.center[1]+h, n.center[2]-h, n.center[2]+h)
end

end # module
