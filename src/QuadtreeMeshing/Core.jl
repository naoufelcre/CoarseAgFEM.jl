module CoreOps

using LinearAlgebra
using ..QuadDefs

export split!, sizing_function_from_distance, logistic_sizing_function, get_leaf_at

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

# Helper to create a sizing function from a generalized distance function
# Returns value for: size <= weight * |distance(x)| + min_h
function sizing_function_from_distance(distance_func::Function, weight::Float64, min_h::Float64)
    return x -> (weight * abs(distance_func(x)) + min_h)
end

# Logistic sizing function for safer coarsening
# h(d) = h_min + (h_max - h_min) / (1 + exp(-k * (|d| - d0)))
function logistic_sizing_function(distance_func::Function, min_h::Float64, max_h::Float64, steepness::Float64, transition_dist::Float64)
    return function(x)
        d = abs(distance_func(x))
        return min_h + (max_h - min_h) / (1.0 + exp(-steepness * (d - transition_dist)))
    end
end

# Forest-aware search
function get_leaf_at(roots::Vector{QuadNode}, point)
    # 1. Identify which root covers the point
    # Optimization: If only 1 root, skip loop
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
        # Point outside domain? Return empty/dummy or error?
        # For robustness in probing, we usually handle bounds checks before calling this.
        # But if called, it implies a bug or float error.
        # Let's return the first root to prevent crash, effectively clamping?
        # Or better: return a sentinel.
        # Existing logic returned `curr` (the leaf or parent where it stopped).
        # Let's return the closest root?
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

# Legacy overload for single root calls (if any remain transiently)
get_leaf_at(root::QuadNode, point) = get_leaf_at([root], point)

end # module
