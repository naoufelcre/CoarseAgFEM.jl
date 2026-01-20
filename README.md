# CoarseAgFEM.jl

Coarse meshs for aggregated finite elements methods of GridapEmbedded.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Naoufel.github.io/CoarseAgFEM.jl/dev)

## API (WIP)

- [Development Documentation](https://Naoufel.github.io/CoarseAgFEM.jl/dev)
    
## Simple Usage

```julia
using Gridap
using CoarseAgFEM

# 1. Define Fine Model & Geometry
model = CartesianDiscreteModel((0,1,0,1), (64,64))
dist(x) = sqrt((x[1]-0.5)^2 + (x[2]-0.5)^2) - 0.25

# 2. Build Coarse Model (One-Liner)
# Returns the Gridap model and the Transfer Operator
coarse_model, op = build_coarse_model(model, dist)

# 3. Visualize
writevtk(coarse_model, "output/coarse_mesh")
```
