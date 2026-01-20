module CoreOps

using LinearAlgebra
using ..QuadDefs

export split!, get_leaf_at

function split!(node::QuadNode, mesh::CoarseMeshBuilder)
    if !is_leaf(node) return end

    h = node.size / 4.0
    new_size = node.size / 2.0
    new_level = node.level + 1

    offsets = [[-h, -h], [ h, -h], [-h,  h], [ h,  h]]

    for i in 1:4
        child = QuadNode(length(mesh.all_nodes)+1, node.center .+ offsets[i], new_size, new_level)
        child.parent = node
        push!(node.children, child)
        push!(mesh.all_nodes, child)
    end
end


# Forest-aware search
function get_leaf_at(roots::Vector{QuadNode}, point)
    #Identify which root covers the point
    #If only 1 root, skip loop
    local root = nothing
    if length(roots) == 1
        root = roots[1]
    else
        for r in roots
            b = get_bounds(r)
            if point[1] >= b[1] && point[1] <= b[2] && point[2] >= b[3] && point[2] <= b[4]
                root = r
                break
            end
        end
    end

    if root === nothing
        return roots[1]
    end

    curr = root
    # Traverse only active nodes logic from original code
    while !is_leaf(curr)
        found = false
        for child in curr.children
            # Basic bound check
            b = get_bounds(child)
            if point[1] >= b[1] && point[1] <= b[2] && point[2] >= b[3] && point[2] <= b[4]
                curr = child
                found = true
                break
            end
        end
        if !found return curr end
    end
    return curr
end

end # module
