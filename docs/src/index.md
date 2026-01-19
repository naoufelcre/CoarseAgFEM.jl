# CoarseAgFEM.jl

Documentation for CoarseAgFEM.jl, a package for robust aggregation finite element methods on Quadtree meshes.

## Installation

```julia
using Pkg
Pkg.add("CoarseAgFEM")
```

## Features

- **CoarseAgFEM**: Robust Aggregation for Finite Element Methods.
- **Quadtree Meshing**: Generate Quadtree meshes from Cartesian grids.
- **Robust Aggregation**: Handle small cut cells with robust aggregation strategies.

## API Reference

```@index
```

```@autodocs
Modules = [CoarseAgFEM]
```

### Exported Functions

- `generate_fine_mesh`
- `bottom_up_coarsening!`
- `classify_leaves!`
- `balance!`
- `pave_mesh`
- `cartesian_to_quadtree`
- `quadtree_to_discrete_model`
- `write_vtk`

### Types

- `QuadMesh`
- `QuadNode`
- `QuadElement`
- `RobustAggregation`
