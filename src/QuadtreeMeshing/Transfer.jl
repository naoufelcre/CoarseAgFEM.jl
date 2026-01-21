module QuadtreeTransferOp

import TransferOperator
using TransferOperator: AbstractTransferOperator
using ..QuadDefs: CoarseMeshBuilder, QuadNode, get_bounds, is_leaf

export QuadtreeTransfer

"""
    QuadtreeTransfer <: AbstractTransferOperator

Transfer operator for Quadtree-based coarse meshes.
"""
struct QuadtreeTransfer <: AbstractTransferOperator
    fine_to_coarse::Array{Int,2}
    n_fine::Tuple{Int,Int}
    n_coarse::Int
    roots::Vector{QuadNode}
    leaf_to_cell::Dict{Int, Vector{Int}}
end

"""
    locate(op::QuadtreeTransfer, x) → coarse_id (actual cell id)

Find coarse cell containing point `x` via tree descent.
Returns the first cell ID found in the leaf.
"""
function TransferOperator.locate(op::QuadtreeTransfer, x)
    for root in op.roots
        leaf_id = _tree_locate_leaf(root, x)
        if leaf_id !== nothing
            # Map leaf_id to its first cell_id
            if haskey(op.leaf_to_cell, leaf_id)
                return op.leaf_to_cell[leaf_id][1]
            end
        end
    end
    return nothing  # Not found
end

function _tree_locate_leaf(node::QuadNode, x)
    b = get_bounds(node)
    if x[1] < b[1] || x[1] > b[2] || x[2] < b[3] || x[2] > b[4]
        return nothing
    end
    if is_leaf(node) && node.is_active
        return node.id
    end
    for child in node.children
        result = _tree_locate_leaf(child, x)
        result !== nothing && return result
    end
    return nothing
end

"""
    restrict(op::QuadtreeTransfer, fine_field) → coarse_field

Restrict fine field to coarse cells (averaging).
Values are averaged per leaf and assigned to all cells in that leaf.
"""
function TransferOperator.restrict(op::QuadtreeTransfer, fine_field::Array{T,2}) where T
    coarse = zeros(T, op.n_coarse)
    counts = zeros(Int, op.n_coarse)
    
    # We first average per leaf
    leaf_sums = Dict{Int, T}()
    leaf_counts = Dict{Int, Int}()
    
    nx, ny = op.n_fine
    for j in 1:ny, i in 1:nx
        lid = op.fine_to_coarse[i,j]
        lid > 0 || continue
        leaf_sums[lid] = get(leaf_sums, lid, zero(T)) + fine_field[i,j]
        leaf_counts[lid] = get(leaf_counts, lid, 0) + 1
    end
    
    # Assign average to all cells in the leaf
    for (lid, s) in leaf_sums
        avg = s / leaf_counts[lid]
        for cell_id in get(op.leaf_to_cell, lid, Int[])
            coarse[cell_id] = avg
        end
    end
    
    return coarse
end

"""
    prolong(op::QuadtreeTransfer, coarse_field) → fine_field

Prolong coarse field to fine grid (piecewise constant by leaf).
"""
function TransferOperator.prolong(op::QuadtreeTransfer, coarse_field::Vector{T}) where T
    nx, ny = op.n_fine
    fill_val = T <: AbstractFloat ? T(NaN) : zero(T)
    fine = fill(fill_val, nx, ny)
    
    for j in 1:ny, i in 1:nx
        lid = op.fine_to_coarse[i,j]
        lid > 0 || continue
        
        if haskey(op.leaf_to_cell, lid)
            # Use the first cell in the leaf (they share the same quadtree node)
            cid = op.leaf_to_cell[lid][1]
            fine[i,j] = coarse_field[cid]
        end
    end
    return fine
end

end # module
