module Visualization

using ..QuadDefs
using Printf

export write_svg, write_vtk

function write_svg(filename, all_elements, root_size=1.0)
    open(filename, "w") do io
        scale = 1000.0
        write(io, "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 $(root_size*scale) $(root_size*scale)' style='background-color:white'>\n")
        write(io, "<style>polygon { stroke: black; stroke-width: $(0.5/scale); vector-effect: non-scaling-stroke; }</style>\n")
        for el in all_elements
            pts_str = join(["$(p[1]*scale), $((root_size - p[2])*scale)" for p in el.nodes], " ") 
            write(io, "<polygon points='$pts_str' fill='$(el.color)' opacity='0.8' />\n")
        end
        write(io, "</svg>")
    end
    println("  Saved SVG to $filename")
end

# ------------------------------------------------------------------------------
# Minimal Legacy VTK Writer (Zero Dependency)
# ------------------------------------------------------------------------------
function write_vtk(filename, elements::Vector{QuadElement}, mesh::Union{QuadMesh, Nothing}=nothing)
    # Gather unique nodes
    unique_nodes = Vector{Float64}[]
    node_map = Dict{Vector{Float64}, Int}()
    
    # Helpers for caching nodes
    function get_node_id(p)
        if !haskey(node_map, p)
            push!(unique_nodes, p)
            node_map[p] = length(unique_nodes) - 1 # 0-based for VTK? 
            # Actually Legacy VTK uses 0-based indexing for connectivity lists
        end
        return node_map[p]
    end

    # Build connectivity
    # VTK uses varying cell types. 
    # Triangle = 5, Quad = 9, Polygon = 7
    
    cells_data = Int[]
    cells_types = Int[]
    
    for el in elements
        n_pts = length(el.nodes)
        push!(cells_data, n_pts)
        for p in el.nodes
            push!(cells_data, get_node_id(p))
        end
        
        if n_pts == 3
            push!(cells_types, 5) # VTK_TRIANGLE
        elseif n_pts == 4
            push!(cells_types, 9) # VTK_QUAD
        else
            push!(cells_types, 7) # VTK_POLYGON
        end
    end
    
    n_nodes = length(unique_nodes)
    n_cells = length(elements)
    list_size = length(cells_data)
    
   open(filename, "w") do io
        write(io, "# vtk DataFile Version 2.0\n")
        write(io, "QuadtreeAggregations Generated Mesh\n")
        write(io, "ASCII\n")
        write(io, "DATASET UNSTRUCTURED_GRID\n")
        
        # POINTS
        write(io, "POINTS $n_nodes double\n")
        for p in unique_nodes
            @printf(io, "%.6f %.6f 0.0\n", p[1], p[2])
        end
        
        # CELLS
        write(io, "\nCELLS $n_cells $list_size\n")
        for el in elements
            ids = [get_node_id(p) for p in el.nodes]
            write(io, "$(length(ids)) " * join(ids, " ") * "\n")
        end
        
        # CELL_TYPES
        write(io, "\nCELL_TYPES $n_cells\n")
        for t in cells_types
            write(io, "$t\n")
        end
        
        # CELL DATA
        write(io, "\nCELL_DATA $n_cells\n")
        
        # 1. Color Type
        write(io, "SCALARS cell_type_color int 1\n")
        write(io, "LOOKUP_TABLE default\n")
        for el in elements
            # Map color string to int for visualization
            val = (el.color == "orange") ? 1 : 0
            write(io, "$val\n")
        end

        # 2. Topology Data (if mesh provided)
        if !isnothing(mesh)
            # Leaf ID
            write(io, "\nSCALARS leaf_id int 1\n")
            write(io, "LOOKUP_TABLE default\n")
            for el in elements
                write(io, "$(el.leaf_id)\n")
            end

            # Level
            write(io, "\nSCALARS level int 1\n")
            write(io, "LOOKUP_TABLE default\n")
            for el in elements
                node = mesh.all_nodes[el.leaf_id]
                write(io, "$(node.level)\n")
            end

            # Parent ID
            write(io, "\nSCALARS parent_id int 1\n")
            write(io, "LOOKUP_TABLE default\n")
            for el in elements
                node = mesh.all_nodes[el.leaf_id]
                pid = isnothing(node.parent) ? 0 : node.parent.id
                write(io, "$pid\n")
            end

            # Is Hanging? (Redundant with color but explicit)
            write(io, "\nSCALARS is_hanging int 1\n")
            write(io, "LOOKUP_TABLE default\n")
            for el in elements
                val = (el.color == "orange") ? 1 : 0
                write(io, "$val\n")
            end
        end
    end
    println("  Saved VTK to $filename")
end

end # module
