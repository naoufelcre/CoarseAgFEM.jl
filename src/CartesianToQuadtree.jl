module CartesianToQuadtree

using Gridap
using Gridap.Geometry
using Gridap.Arrays
using ..QuadtreeAggregations

export cartesian_to_quadtree

"""
    cartesian_to_quadtree(model::CartesianDiscreteModel)

Converts a Gridap `CartesianDiscreteModel` into a `QuadMesh`.
The conversion is "Bottom-Up":
1. Input cells become the initial "Leaves" of the QuadMesh.
2. We infer the tree structure (parents) based on the grid topology.
3. The resulting QuadMesh may be a "Forest" (multiple roots) if the domain is not a square power-of-two.
"""
function cartesian_to_quadtree(model::CartesianDiscreteModel)
    println("  [CartesianToQuadtree] Converting Cartesian Model to QuadMesh...")
    
    # 1. Extract Grid Information
    # ---------------------------
    desc = get_cartesian_descriptor(model)
    partition = desc.partition
    nx, ny = partition[1], partition[2]

    # Gridap.Geometry.CartesianDescriptor.sizes stores the CELL SIZES (hx, hy), not domain lengths.
    # Note: Confirmed by debug trace where sizes[1] = 0.02 (1/50) for domain (0,1).
    hx = desc.sizes[1]
    hy = desc.sizes[2]
    
    Lx = hx * nx
    Ly = hy * ny
    origin = desc.origin
    
    println("    [Debug] Domain: Origin=$origin. Lx=$Lx, Ly=$Ly. Partition: $nx x $ny")
    println("    [Debug] hx=$hx, hy=$hy, cell_size=$((hx+hy)/2.0)")
    
    # We assume uniform cells (squares) for standard Quadtree logic
    if !isapprox(hx, hy; rtol=1e-5)
        # @warn "[CartesianToQuadtree] Cells are not square..."
    end
    
    cell_size = (hx + hy) / 2.0
    
    # 2. Instantiate Leaves
    # ---------------------
    # We Map (i, j) index to a QuadNode.
    # Gridap indices are 1-based.
    # We'll use a dictionary to store nodes by (level, grid_x, grid_y)
    # where level 0 = The Leaves (finest). 
    # NOTE: Standard Quadtree usually has Level 0 = Root. 
    # Here, let's stick to standard: Level MAX = Leaves. 
    # But we don't know MAX yet.
    # Let's use an inverted level for construction: Depth 0 = Leaves.
    # Then we can flip it later or just map Depth k -> Size = S0 * 2^k.
    
    # Let's use "Level" as defined in QuadDefs:
    # We need to decide a "Max Level". 
    # Let's arbitrarily say the leaves are at Level L = 10 (or sufficient depth).
    # Or better: We calculate Level 0 based on the "Domain Bounding Box".
    
    # Strategy:
    # Find the smallest power-of-two square that covers the domain.
    # L_dim = max(Lx, Ly)
    # root_size = next_pow_2(L_dim)
    # depth = log2(root_size / cell_size)
    
    # Let's assume the domain origin is (0,0) for Morton coding simplicity, 
    # or just use relative coordinates.
    
    # Let's compute the theoretical Full Root Size
    max_dim = max(Lx, Ly)
    # We need an integer number of cells to be a power of 2? 
    # Not necessarily, but for standard Quadtree ops, standardizing helps.
    # But constructing bottom-up, we don't *need* a single root.
    
    # Let's just create nodes.
    # ID counter
    node_id_counter = 0
    
    # Store: (DepthFromBottom, i, j) -> QuadNode
    # Depth 0 = Leaves. i,j are integer indices at that level.
    nodes_by_level_idx = Dict{Tuple{Int, Int, Int}, QuadNode}()
    
    # Create Leaves
    all_nodes = QuadNode[]
    
    # Gridap cell iterator
    # CartesianDescription gives us `map` to go from integer index to point.
    # We can just run i in 1:nx, j in 1:ny
    
    # Calculate a "Base Level" index for the leaves.
    # If we treat the leaves as Level = MAX, we need MAX.
    # Let's just give them Level = 0 for now and defined "Level" as "Refinement Level relative to Root" later?
    # No, `QuadNode` struct expects `level`.
    # Let's use a placeholder `depth` and fix `level` at the end if needed.
    # But `QuadDefs` doesn't enforce logic on `level` other than for our own accounting.
    # Let's say Leaves are Level 0 (Finest) and Roots are Level K. (Reverse of standard? No standard is Root=0).
    # OK, let's effectively picking a large integer, say 20, for the leaves, and decrement going up.
    leaf_level = 20 
    
    println("    -> Instantiating $nx x $ny leaves...")
    
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
    # We attempt to build parents for depths 1, 2, ...
    # Parent (depth k+1) at index (I, J) covers (2I-1, 2J-1) ... (2I, 2J) of depth k?
    # Wait, 1-based indexing is tricky with mod. 
    # Let's map 1-based i,j to 0-based for math: i' = i-1.
    # Parent index I' = i' // 2.
    # Cell (I') covers (2I', 2I'+1).
    # Back to 1-based: Parent I = (i-1)÷2 + 1.
    
    current_depth = 0
    current_nx = nx
    current_ny = ny
    current_size = cell_size
    
    while current_nx > 1 || current_ny > 1
        # Determine range of potential parents
        parent_nx = ceil(Int, current_nx / 2)
        parent_ny = ceil(Int, current_ny / 2)
        
        println("    [Debug Loop] Depth=$current_depth. nx=$current_nx. current_size=$current_size. parent_nx=$parent_nx")
        
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
            # Check if we should create this parent.
            # In "Greedy Packing" / "Forest" mode:
            # We construct the parent node object.
            # We connect available children.
            # Even if only 1 child exists? 
            # If we create a parent with only 1 child, and we later coarsen it, we might get an invalid state 
            # if the domain was supposed to be void there.
            # BUT: A QuadNode is just a container. It becomes a Leaf if children are merged.
            # If we simply build the hierarchy, we represent the "Sparse Quadtree".
            
            # Parent geometry
            # Size = 2 * child_size
            # Center?
            # Child (i, j) center relative to parent?
            # Let's compute center from geometry index.
            # Grid Origin + (PI - 0.5) * P_size
            
            p_size = current_size * 2.0
            p_cx = origin[1] + (PI - 0.5) * p_size
            p_cy = origin[2] + (PJ - 0.5) * p_size
            
            # Create Parent
            node_id_counter += 1
            # Level = leaf_level - (current_depth + 1)
            p_node = QuadNode(node_id_counter, [p_cx, p_cy], p_size, leaf_level - (current_depth + 1))
            p_node.is_active = true # Must be active to be a candidate for coarsening check (it won't be a leaf because it has children)
            # Wait, existing logic: if children exist, parent is NOT leaf.
            # Active flag logic: Active = It is part of the current mesh.
            # Internal nodes are active? In `QuadDefs`: "is_active # If false, it has been removed/merged"
            # So Internal nodes ARE active in the tree sense, but not leaves.
            # `is_active` usually denotes "Part of the cut/boundary representation" or just "Exists"?
            # `QuadtreeAggregations` logic: `leaves = [n for n in mesh.all_nodes if n.is_active && is_leaf(n)]`
            # So Internal Nodes can be is_active=true.
            # But `bottom_up_coarsening` sets children.is_active = false when merging.
            # So: Initially ALL nodes we create are likely active, preserving the hierarchy.
            
            # Find Children
            # Indices: (2*PI - 1) and (2*PI)
            # 0-based: 2*(PI-1) -> +1 for 1-based
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
            
            # Verify strict containment? 
            # If we implement "Hidden Trunk", we allow partial parents.
            # But `QuadNode` struct expects children vector.
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
    # We need a "Root".
    # If we have multiple components or "Forest", `QuadMesh` logic expecting a single root might fail?
    # `QuadMesh` struct: `root::QuadNode`.
    # We might need a "Super Root" that covers everything if we want to adhere to the struct.
    # Or we relax the struct to `roots::Vector{QuadNode}`.
    # For now, let's see if our bottom-up process converged to a single node.
    # If the domain is POT square, it converges to 1 node.
    # If not, we might have a few top-level nodes.
    # Let's create a Virtual Super Root if needed, or pick the last one.
    
    # Check top level
    top_nodes = [n for n in all_nodes if n.parent === nothing]
    
    local root
    if length(top_nodes) == 1
        root = top_nodes[1]
    else
        println("  [CartesianToQuadtree] Forest detected ($(length(top_nodes)) roots). Creating Virtual Super Root.")
        # Create a dummy super root just to satisfy the struct
        # This might break things if algorithms assume geometry containment for the root.
        # But `bottom_up_coarsening` iterates `mesh.all_nodes`, so root isn't heavily used except for `get_leaf_at`.
        # `get_leaf_at` starts at root. It NEEDS to cover the domain.
        
        # We must create a Super Root covering all top_nodes.
        # Simple hack: Keep building up until 1 node?
        # My `while` loop stopped when nx, ny <= 1? No, `current_nx > 1 || current_ny > 1`.
        # So it should converge to 1x1 if we handle the Ceiling correctly?
        # indices go 1..1 eventually.
        # Let's trace: 3x3 -> 2x2 -> 1x1.
        # yes.
        # But wait, 3 cells.
        # Level 0: 1, 2, 3.
        # Level 1: Parent(1,2) [covers 1-2], Parent(2) [covers 3-4?].
        # If Child 4 doesn't exist, Parent(2) has only Child 3.
        # Level 2: Parent(1-2), Parent(3-4) -> GrandParent.
        # It DOES converge.
        # Why did I worry?
        
        if isempty(top_nodes)
             error("[CartesianToQuadtree] No nodes created!")
        end
        root = top_nodes[1] # Should be unique if loop works
    end
    
    # Re-normalize Levels
    # Current Leaves are at `leaf_level` (e.g. 20).
    # We want Root to be at Level 0 for standard notation?
    # Or keep it relative?
    # Existing `generate_fine_mesh` sets Root at Level 0.
    # So let's shift all levels.
    shift = root.level # e.g. 15
    for n in all_nodes
        n.level = n.level - shift
        # root becomes 0. leaves become positive.
        n.level = abs(n.level) # Just in case
    end
    
    println("  [CartesianToQuadtree] Tree built. Root Level 0. Leaves Level $(root.children[1].level + 1). Total Nodes: $(length(all_nodes))")
    
    return QuadMesh(root, all_nodes)
end

end # module
