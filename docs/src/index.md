# CoarseAgFEM.jl

Coarse meshs for aggregated finite elements methods of GridapEmbedded.jl

## API (WIP)

## API Reference

```@index
```

```@autodocs
Modules = [
    CoarseAgFEM, 
    CoarseAgFEM.CartesianToQuadtree, 
    CoarseAgFEM.RobustAgFEM, 
    CoarseAgFEM.GridapIntegration,
    CoarseAgFEM.QuadtreeAggregations,
    CoarseAgFEM.QuadtreeAggregations.Coarsening
]
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
