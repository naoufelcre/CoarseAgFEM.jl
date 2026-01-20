module CoarseAgFEM

using Gridap
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Helpers
using Gridap.ReferenceFEs
using GridapEmbedded
using GridapEmbedded.CSG
using GridapEmbedded.Interfaces

# Re-export key types
export QuadMesh, QuadNode, QuadElement, NodeStatus, INTERIOR, EXTERIOR, CUT, UNDEFINED
export generate_fine_mesh, bottom_up_coarsening!, classify_leaves!, balance!, pave_mesh, cartesian_to_quadtree
export quadtree_to_discrete_model, build_coarse_model
# export sizing_function... (deprecated)
export RobustAggregation
export write_vtk
export TransferOperator

include("QuadtreeMeshing/QuadtreeMeshing.jl")
include("CartesianToQuadtree.jl")
include("GridapIntegration.jl")
include("RobustAgFEM.jl")
include("TransferOperator.jl")

# Re-export from submodules
using .QuadtreeMeshing
using .CartesianToQuadtree
using .GridapIntegration
using .RobustAgFEM
using .TransferOperators

"""
    build_coarse_model(fine_model::CartesianDiscreteModel, level_set_func::Function; parameters...)

High-level API to construct a CoarseAgFEM model.
Returns `(coarse_model, transfer_operator)`.
"""
function build_coarse_model(fine_model::CartesianDiscreteModel, level_set_func::Function; 
                            max_coarsening_factor=8, 
                            buffer_width=0.0)
    
    # 1. Instantiate Builder (Internal)
    builder = cartesian_to_quadtree(fine_model)
    
    # 2. Refine/Coarsen (Internal logic)
    classify_leaves!(builder, level_set_func; buffer_width=buffer_width)
    bottom_up_coarsening!(builder; max_coarsening_factor=max_coarsening_factor)
    balance!(builder)
    
    # 3. Finalize
    # Extract Elements (Paving)
    elements = pave_mesh(builder)
    
    # Convert to Gridap Model
    coarse_model, leaf_to_cell_ids = quadtree_to_discrete_model(elements)
    
    # Compute Transfer Operator
    # We need to map fine cells to the final coarse cells.
    # We can rely on generate_fine_to_coarse_map -> gives us Leaf ID.
    # Then leaf_to_cell_ids -> gives us Gridap Cell IDs.
    
    # Is fine_model uniform?
    desc = get_cartesian_descriptor(fine_model)
    nx = desc.partition[1]
    ny = desc.partition[2]
    
    # Generate Map matching the exact partition
    raw_map = generate_fine_to_coarse_map(builder, (nx, ny))
    
    op = TransferOperator(raw_map, (nx, ny), length(elements)) # n_coarse approximation
    
    return coarse_model, op
end

end # module
