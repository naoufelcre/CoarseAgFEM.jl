module TransferOperators

export TransferOperator

"""
    struct TransferOperator

Represents the transfer operator (prolongation/restriction) between a fine grid and a coarse mesh.
Currently stores the aggregation map (piecewise constant basis).

Fields:
- `fine_to_coarse::Array{Int, 2}`: Map where `fine_to_coarse[i,j]` is the ID of the coarse cell containing fine cell `(i,j)`.
- `n_fine::Tuple{Int, Int}`: Dimensions of the fine grid.
- `n_coarse::Int`: Total number of coarse cells.
"""
struct TransferOperator
    fine_to_coarse::Array{Int, 2}
    n_fine::Tuple{Int, Int}
    n_coarse::Int
end

# TODO: Add methods for restrict/prolong

end
