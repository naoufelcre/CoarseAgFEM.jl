# CoarseAgFEM.jl

Coarse meshs for aggregated finite elements methods of GridapEmbedded.jl

## API (WIP)

## API Reference

```@index
```

```@autodocs
Modules = [
    CoarseAgFEM, 
    CoarseAgFEM.MeshCoarsening, 
    CoarseAgFEM.RobustAgFEM, 
    CoarseAgFEM.GridapIntegration,
    CoarseAgFEM.QuadtreeMeshing,
    CoarseAgFEM.QuadtreeMeshing.QuadDefs,
    CoarseAgFEM.QuadtreeMeshing.Coarsening,
    CoarseAgFEM.TransferOperators
]
```

### Exported Functions

- `generate_fine_mesh`
- `build_coarse_model`
- `bottom_up_coarsening!`
- `classify_leaves!`
- `balance!`
- `pave_mesh`
- `initialize_builder`
- `quadtree_to_discrete_model`
- `write_vtk`

### Types

- `CoarseMeshBuilder`
- `QuadNode`
- `QuadElement`
- `TransferOperator`
- `RobustAggregation`
