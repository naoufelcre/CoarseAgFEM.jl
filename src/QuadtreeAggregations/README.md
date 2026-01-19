# QuadtreeAggregations Plugin

A standalone, zero-dependency Julia module for generating adaptive quadtree meshes via bottom-up aggregation.

## Usage Overview

The workflow consists of four standard phases:
1.  **Generation**: Create a uniform fine grid.
2.  **Coarsening**: Aggregates cells based on a customizable sizing function.
3.  **Balancing**: Enforces 2:1 sizing constraint between neighbors.
4.  **Paving**: Resolves hanging nodes using template triangulation.

## API Reference

### 1. Mesh Generation
```julia
mesh = generate_fine_mesh(L_size, max_level)
```

### 2. Coarsening (Aggregation) with Custom Geometry

The coarsening phase now accepts a general `sizing_function` that defines the target element size at any given point in the domain.

```julia
# 1. Define your geometry (e.g., signed distance function)
dist_func(x) = norm(x .- center) - R

# 2. Create a sizing function helper
#    Returns a function: x -> weight * |dist(x)| + min_h
sizing = sizing_function_from_distance(dist_func, weight, min_h)

# 3. Run coarsening
bottom_up_coarsening!(mesh, max_level, sizing)
```
-   `sizing_function`: A function `f(x::Vector) -> Float64` that returns the MAXIMUM allowed element size at point `x`.
-   `sizing_function_from_distance`: A helper to easily create this from a distance field.

### 3. Balancing
```julia
balance!(mesh)
```

### 4. Paving (Meshing)
```julia
elements = pave_mesh(mesh)
```

### 5. Export
```julia
write_svg(filename, elements)
write_vtk(filename, elements)
```

## Extensibility

### Custom Sizing Laws
You are not limited to distance fields. You can provide ANY function `f(x)` for `bottom_up_coarsening!`. For example, you could adapt based on an error estimate, a density field, or a complex CSG geometry.

```julia
my_sizing(x) = (x[1] < 0.5) ? 0.01 : 0.1
bottom_up_coarsening!(mesh, 6, my_sizing)
```
