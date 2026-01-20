module Coarsening

using ..QuadDefs
using ..CoreOps
using Printf

export generate_fine_mesh, bottom_up_coarsening!, classify_leaves!, generate_fine_to_coarse_map

# Legacy function kept for compatibility if needed, but not used in new flow
function generate_fine_mesh(L, max_level)
    # ... (Keep minimal or error out? Let's just keep it simple or empty)
    error("generate_fine_mesh is deprecated. Use cartesian_to_quadtree.")
end

"""
    classify_leaves!(mesh::CoarseMeshBuilder, level_set_func::Function; buffer_width::Float64=0.0)

Classifies all active leaves as INTERIOR, EXTERIOR, CUT, or BUFFER.
BUFFER status is assigned if the cell is not CUT coverage, but is within `buffer_width` of the interface.
"""
function classify_leaves!(mesh::CoarseMeshBuilder, level_set_func::Function; buffer_width::Float64=0.0)
    println("  Phase 0: Classifying Leaves (Buffer Width=$buffer_width)...")
    count_int = 0
    count_ext = 0
    count_cut = 0
    count_buf = 0
    
    leaves = [n for n in mesh.all_nodes if n.is_active && is_leaf(n)]
    
    for leaf in leaves
        # Evaluate phi at 4 comers
        b = get_bounds(leaf) 
        corners = [
            [b[1], b[3]], # SW
            [b[2], b[3]], # SE
            [b[2], b[4]], # NE
            [b[1], b[4]]  # NW
        ]
        
        vals = [level_set_func(p) for p in corners]
        
        has_pos = any(v -> v > 0, vals)
        has_neg = any(v -> v < 0, vals)
        
        # Min absolute value represents approx distance to interface
        min_abs_dist = minimum(abs.(vals))
        
        if has_pos && has_neg
            leaf.status = CUT
            count_cut += 1
        elseif min_abs_dist < buffer_width
            # Close to interface, but strictly on one side
            leaf.status = BUFFER
            count_buf += 1
        elseif has_neg # Interior and Far
            leaf.status = INTERIOR
            count_int += 1
        else # Exterior and Far
            leaf.status = EXTERIOR
            count_ext += 1
        end
    end
    println("    -> Classification Results: Interior=$count_int, Exterior=$count_ext, Cut=$count_cut, Buffer=$count_buf")
end

"""
    bottom_up_coarsening!(mesh::CoarseMeshBuilder; max_coarsening_factor::Int=8)

Greedy packing: Merges 4 siblings if they are HOMOGENEOUS (all Interior or all Exterior).
Cut cells and Mixed cells are never merged.
Coarsening stops if the parent cell size would exceed `max_coarsening_factor * h_min` (e.g. 8*h).
"""
function bottom_up_coarsening!(mesh::CoarseMeshBuilder; max_coarsening_factor::Int=8)
    println("  Phase 1: Bottom-Up Homogeneity Coarsening (Max Factor = $max_coarsening_factor)...")
    
    # Determine max level (finest)
    max_level = maximum(n.level for n in mesh.all_nodes)
    
    # Calculate min allowed level
    # 2^k = factor => k = log2(factor)
    # parent_level = max_level - k
    # If factor=8, k=3. We allow parents at level max-3.
    # We DO NOT allow parents at level max-4.
    
    depth_limit = Int(floor(log2(max_coarsening_factor)))
    min_allowed_level = max_level - depth_limit
    
    println("    -> Max Level: $max_level. Min Allowed Level: $min_allowed_level (Limit for Size <= $(max_coarsening_factor)h).")
    
    merges_count = 0
    
    # Iterate levels Bottom-Up
    for l in (max_level-1):-1:0
        
        # STOP if we are trying to form parents at a forbidden level
        if l < min_allowed_level
            # We are at level l. Candidates are parents at level l.
            # If l < min_allowed, these parents are too big.
            continue 
        end
        
        candidates = [n for n in mesh.all_nodes if n.is_active && n.level == l]
        
        # We process logical parents.
        # Note: 'candidates' here are the parents we want to potentially form.
        # But wait, in the 'Forest' construction, parents are ALREADY instantiated and is_active=true.
        # So we check if *their children* are mergeable.
        
        for parent in candidates
            if is_leaf(parent)
                continue
            end
            
            # 1. Must have all 4 children?
            # In optimal packing, we only merge full blocks?
            # Or do we handle boundary robustness?
            # STRICT RULE: Must have 4 children to form a valid larger square.
            if length(parent.children) != 4
                # Cannot merge if not a full quad
                continue
            end
            
            # 2. Check Children Status
            # Children must be Active Leaves (or previously merged nodes which are now leaves)
            children_ready = true
            first_status = parent.children[1].status
            
            for child in parent.children
                if !child.is_active || !is_leaf(child) # Must be a current leaf
                    children_ready = false
                    break
                end
                
                # Homogeneity Check
                if child.status != first_status
                    children_ready = false
                    break
                end
                
                # Check allowed types
                if child.status == CUT || child.status == BUFFER || child.status == UNDEFINED
                    children_ready = false # Never merge cuts or buffer cells
                    break
                end
            end
            
            if children_ready
                # MERGE!
                # 1. Mark children inactive
                for child in parent.children
                    child.is_active = false
                end
                # 2. Inherit status
                parent.status = first_status
                # 3. Disconnect children to make parent a leaf
                empty!(parent.children)
                
                merges_count += 1
            else
                # If we don't merge, the parent status remains UNDEFINED (virtual container)
                # unless it was already set? No, initialized UNDEFINED.
                # It just acts as a branch.
            end
        end
    end
    println("  -> Merged $merges_count parents.")
end

"""
    generate_fine_to_coarse_map(mesh::CoarseMeshBuilder, max_level::Int)

Generates a mapping between the uniform fine grid (at `max_level`) and the current
coarse QuadMesh.

Returns:
    - fine_map: Matrix{Int} of size (2^max_level, 2^max_level). 
      fine_map[i, j] contains the ID of the Leaf Node that covers fine cell (i, j).
    
This serves as a mapping operator M: Ω_h -> Ω_H.
"""
# Legacy/Convenience: If only max_level is provided, assume POT
function generate_fine_to_coarse_map(mesh::CoarseMeshBuilder, max_level::Int)
    N = 2^max_level
    return generate_fine_to_coarse_map(mesh, (N, N))
end

"""
    generate_fine_to_coarse_map(mesh::CoarseMeshBuilder, partition::Tuple{Int, Int})

Generates a mapping between the uniform fine grid (defined by `partition`) and the current
coarse CoarseMeshBuilder.
"""
function generate_fine_to_coarse_map(mesh::CoarseMeshBuilder, partition::Tuple{Int, Int})
    nx, ny = partition
    println("  Phase 2: Generating Fine-to-Coarse Mapping (Grid: $nx x $ny)...")
    
    # We assume domain is (0, Lx) x (0, Ly) ??
    # Actually, CoarseMeshBuilder doesn't store explicit global Domain Bounds easily accessible?
    # Leaf 1 center and size gives us info. 
    # But usually Quadtree assumes relative coords or handled during construction.
    # CartesianToQuadtree uses `desc.origin` and `desc.sizes`.
    
    # We need to know hx, hy.
    # We assume 'mesh' covers the same domain as the partition.
    # But we don't have the domain size here! 
    # `QuadNode` has `center` and `size`. 
    # Bounding box of all nodes?
    
    # Let's compute global bounds from roots
    # Assumes roots cover the domain.
    # Or just iterate all leaves and compute BB.
    
    xmin = Inf; xmax = -Inf; ymin = Inf; ymax = -Inf
    active_leaves = [n for n in mesh.all_nodes if n.is_active && is_leaf(n)]
    
    for leaf in active_leaves
        b = get_bounds(leaf)
        xmin = min(xmin, b[1])
        xmax = max(xmax, b[2])
        ymin = min(ymin, b[3])
        ymax = max(ymax, b[4])
    end
    
    Lx = xmax - xmin
    Ly = ymax - ymin
    
    hx = Lx / nx
    hy = Ly / ny
    
    # Resulting map
    fine_map = zeros(Int, nx, ny)
    
    count_filled = 0
    eps = 1e-10 * min(hx, hy)
    
    for leaf in active_leaves
        b = get_bounds(leaf)
        
        # Calculate index ranges relative to (xmin, ymin)
        # 1-based indices
        
        ix_start = floor(Int, (b[1] - xmin) / hx + eps) + 1
        ix_end   = floor(Int, (b[2] - xmin) / hx - eps) + 1
        
        iy_start = floor(Int, (b[3] - ymin) / hy + eps) + 1
        iy_end   = floor(Int, (b[4] - ymin) / hy - eps) + 1
        
        # Clamp
        ix_start = clamp(ix_start, 1, nx)
        ix_end   = clamp(ix_end,   1, nx)
        iy_start = clamp(iy_start, 1, ny)
        iy_end   = clamp(iy_end,   1, ny)
        
        for y in iy_start:iy_end
            for x in ix_start:ix_end
                fine_map[x, y] = leaf.id
                count_filled += 1
            end
        end
    end
    
    coverage_pct = (count_filled / (nx*ny)) * 100
    println("    -> Mapping generated. Covered $(count_filled)/$(nx*ny) fine cells ($(Printf.@sprintf("%.2f", coverage_pct))%).")
    
    return fine_map
end

end # module
