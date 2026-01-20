# CoarseAgFEM.jl

Coarse meshs for aggregated finite elements methods of GridapEmbedded.jl

[![API](https://img.shields.io/badge/docs-dev-blue.svg)](https://Naoufel.github.io/CoarseAgFEM.jl)

## API (WIP)

- [Development Documentation](https://Naoufel.github.io/CoarseAgFEM.jl/dev)

## Simple Usage

```julia
using Gridap
using CoarseAgFEM

# Create a uniform grid and level set function
model = CartesianDiscreteModel((0,1,0,1), (64,64))
level_set(x) = sqrt((x[1]-0.5)^2 + (x[2]-0.5)^2) - 0.25

# Build Coarse Model returns the Gridap model and the Transfer Operator
coarse_model, op = build_coarse_model(model, dist)

# Ready to go !

#Visualize
writevtk(coarse_model, "output/coarse_mesh")
```
