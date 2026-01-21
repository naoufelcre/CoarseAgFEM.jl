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
export CoarseMeshBuilder, QuadNode, QuadElement, NodeStatus, INTERIOR, EXTERIOR, CUT, BUFFER, UNDEFINED
export generate_fine_mesh, bottom_up_coarsening!, classify_leaves!, balance!, pave_mesh, initialize_builder, get_bounds
export quadtree_to_discrete_model, build_coarse_model
# export sizing_function... (deprecated)
export RobustAggregation
export write_vtk, write_svg
export QuadtreeTransfer

include("QuadtreeMeshing/QuadtreeMeshing.jl")
include("GridapIntegration.jl")
include("MeshCoarsening.jl")
include("RobustAgFEM.jl")

# Re-export from submodules
using .QuadtreeMeshing
using .MeshCoarsening
using .GridapIntegration
using .RobustAgFEM

# Resolve potential namespace conflicts (e.g. CUT with GridapEmbedded)
const INTERIOR  = QuadtreeMeshing.INTERIOR
const EXTERIOR  = QuadtreeMeshing.EXTERIOR
const CUT       = QuadtreeMeshing.CUT
const BUFFER    = QuadtreeMeshing.BUFFER
const UNDEFINED = QuadtreeMeshing.UNDEFINED

end # module
