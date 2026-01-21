module QuadtreeMeshing

# Sub-modules (pointing to subdirectory)
include("QuadDefs.jl")
include("Core.jl")
include("Coarsening.jl")
include("Balancing.jl")
include("Paving.jl")
include("Visualization.jl")
include("Transfer.jl")

using .QuadDefs
using .CoreOps
using .Coarsening
using .Balancing
using .Paving
using .Visualization
using .QuadtreeTransferOp

# Re-export public API
export CoarseMeshBuilder, QuadNode, QuadElement, get_leaf_at, get_bounds, NodeStatus, INTERIOR, EXTERIOR, CUT, BUFFER, UNDEFINED
export generate_fine_to_coarse_map
export bottom_up_coarsening!, classify_leaves!
export balance!
export pave_mesh
export write_svg, write_vtk
export QuadtreeTransfer

end # module
