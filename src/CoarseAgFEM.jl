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
export quadtree_to_discrete_model
# export sizing_function... (deprecated)
export RobustAggregation
export write_vtk

include("QuadtreeAggregations/QuadtreeAggregations.jl")
include("CartesianToQuadtree.jl")
include("GridapIntegration.jl")
include("RobustAgFEM.jl")

# Re-export from submodules
using .QuadtreeAggregations
using .CartesianToQuadtree
using .GridapIntegration
using .RobustAgFEM

end # module
