module CartesianToQuadtree

using Gridap
using Gridap.Geometry
using Gridap.Arrays
using ..QuadtreeMeshing

export cartesian_to_quadtree

"""
    cartesian_to_quadtree(model::CartesianDiscreteModel)

Converts a Gridap `CartesianDiscreteModel` into a `CoarseMeshBuilder`.
The conversion is "Bottom-Up":
1. Input cells become the initial "Leaves" of the QuadMesh.
2. We infer the tree structure (parents) based on the grid topology.
3. The resulting CoarseMeshBuilder may be a "Forest" (multiple roots) if the domain is not a square power-of-two.
"""
function cartesian_to_quadtree(model::CartesianDiscreteModel)

    
    # 1. Extract Grid Information
    # ---------------------------
    desc = get_cartesian_descriptor(model)
    partition = desc.partition
    nx, ny = partition[1], partition[2]

    # Gridap.Geometry.CartesianDescriptor.sizes stores the CELL SIZES (hx, hy)
    hx = desc.sizes[1]
    hy = desc.sizes[2]
    
    Lx = hx * nx
    Ly = hy * ny
    origin = desc.origin
    
    # We assume uniform cells (squares) for standard Quadtree logic
    if !isapprox(hx, hy; rtol=1e-5)
        # Non-square cells are not strictly supported by standard Quadtree logic,
        # but may work if the aspect ratio is handled elsewhere.
    end
    
    cell_size = (hx + hy) / 2.0
    
    # 2. Instantiate Leaves
    # ---------------------
    # We instantiate leaves at an arbitrary high level (depth) and later normalize
    # so that the root ends up at Level 0.
    node_id_counter = 0
    
    # Store: (DepthFromBottom, i, j) -> QuadNode
    # Depth 0 = Leaves. i,j are integer indices at that level.
    nodes_by_level_idx = Dict{Tuple{Int, Int, Int}, QuadNode}()
    
    all_nodes = QuadNode[]
    
    leaf_level = 20 # Arbitrary high number, will be shifted later 
    

    
    for j in 1:ny
        for i in 1:nx
            # Center computation
            # x = origin[1] + (i - 0.5)*hx
            cx = origin[1] + (i - 0.5)*hx
            cy = origin[2] + (j - 0.5)*hy
            center = [cx, cy]
            
            node_id_counter += 1
            node = QuadNode(node_id_counter, center, cell_size, leaf_level)
            
            # Store in map
            nodes_by_level_idx[(0, i, j)] = node
            push!(all_nodes, node)
        end
    end
    
    # 3. Bottom-Up Construction
    # -------------------------
    # Recursively group 2x2 blocks of nodes into parent nodes until the entire domain is covered.
    
    current_depth = 0
    current_nx = nx
    current_ny = ny
    current_size = cell_size
    
    while current_nx > 1 || current_ny > 1
        # Determine range of potential parents
        parent_nx = ceil(Int, current_nx / 2)
        parent_ny = ceil(Int, current_ny / 2)
        

        
        # We prefer to use a Set of parents to create, based on existing children
        potential_parents = Set{Tuple{Int, Int}}()
        
        # Iterate existing nodes at current_depth to flag parents
        for ((d, i, j), node) in nodes_by_level_idx
            if d != current_depth continue end
            
            parent_i = div(i - 1, 2) + 1
            parent_j = div(j - 1, 2) + 1
            push!(potential_parents, (parent_i, parent_j))
        end
        
        if isempty(potential_parents) break end
        
        created_parents = 0
        
        for (PI, PJ) in potential_parents
            # Construct Parent Node
            # The parent covers a 2x2 area in the child grid.
            
            p_size = current_size * 2.0
            p_cx = origin[1] + (PI - 0.5) * p_size
            p_cy = origin[2] + (PJ - 0.5) * p_size
            
            # Create Parent
            node_id_counter += 1
            p_node = QuadNode(node_id_counter, [p_cx, p_cy], p_size, leaf_level - (current_depth + 1))
            p_node.is_active = true 
            
            # Find Children
            # Calculate the range of child indices that belong to this parent
            c_i_start = (PI - 1) * 2 + 1
            c_j_start = (PJ - 1) * 2 + 1
            
            children_found = QuadNode[]
            
            for d_j in 0:1
                for d_i in 0:1
                    ci = c_i_start + d_i
                    cj = c_j_start + d_j
                    key = (current_depth, ci, cj)
                    if haskey(nodes_by_level_idx, key)
                        child = nodes_by_level_idx[key]
                        child.parent = p_node
                        push!(children_found, child)
                    end
                end
            end
            
            # Assign children
            p_node.children = children_found
            
            # Store Parent and Add to Pool
            nodes_by_level_idx[(current_depth + 1, PI, PJ)] = p_node
            push!(all_nodes, p_node)
            created_parents += 1
        end
        
        # Advance
        current_depth += 1
        current_nx = parent_nx
        current_ny = parent_ny
        current_size *= 2.0
    end
    
    # 4. Construct QuadMesh
    # ---------------------
    # The bottom-up process will produce one or more roots (a Forest) depending on the
    # domain dimensions (POT vs non-POT).
    
    roots = [n for n in all_nodes if n.parent === nothing]
    
    if isempty(roots)
         error("[CartesianToQuadtree] No nodes created!")
    end

    # Re-normalize Levels
    # Shift levels so that the Roots are at Level 0.
    # Note: If it's a forest, all roots "should" be at the same conceptual top level 
    # if constructed from a uniform grid, but technically they might vary if the domain is irregular.
    # For simplicity, we just take the max level among roots as the baseline.
    
    # Actually, in our bottom-up build, we assigned `leaf_level - depth`.
    # The 'depth' reached might differ if some parts of the forest are deeper?
    # No, we iterate uniformly. The roots will be at: leaf_level - max_depth_reached.
    # So `roots[1].level` is fine.
    
    shift = roots[1].level
    for n in all_nodes
        n.level = abs(n.level - shift) 
    end
    
    return CoarseMeshBuilder(roots, all_nodes)
end

end # module
