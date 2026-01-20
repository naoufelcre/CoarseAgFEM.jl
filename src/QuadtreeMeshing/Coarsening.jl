module Coarsening

using ..QuadDefs
using ..CoreOps
using Printf

export bottom_up_coarsening!, classify_leaves!, generate_fine_to_coarse_map

"""
    classify_leaves!(mesh::CoarseMeshBuilder, level_set_func::Function; buffer_width::Float64=0.0)

Classifies all active leaves as INTERIOR, EXTERIOR, CUT, or BUFFER.
BUFFER status is assigned if the cell is not CUT coverage, but is within `buffer_width` of the interface.
"""
function classify_leaves!(mesh::CoarseMeshBuilder, level_set_func::Function; buffer_width::Float64=0.0)
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
Remark, c'est peut-être mieux de donner l'exposant de la puissance de deux direct.
"""
function bottom_up_coarsening!(mesh::CoarseMeshBuilder; max_coarsening_factor::Int=8)

    # Determine max level (finest)
    max_level = maximum(n.level for n in mesh.all_nodes)

    depth_limit = Int(floor(log2(max_coarsening_factor)))
    min_allowed_level = max_level - depth_limit

    println("    -> Max Level: $max_level. Min Allowed Level: $min_allowed_level (Limit for Size <= $(max_coarsening_factor)h).")

    merges_count = 0

    # Iterate levels Bottom-Up
    for l in (max_level-1):-1:0

        # STOP if we are trying to form parents at a forbidden level
        if l < min_allowed_level
            continue
        end

        candidates = [n for n in mesh.all_nodes if n.is_active && n.level == l]


        # we check if childrens are mergeable for each parent.

        for parent in candidates
            if is_leaf(parent)
                continue
            end

            # STRICT RULE: Must have 4 children to form a valid larger square.
            if length(parent.children) != 4
                continue
            end

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
                #Mark inactive
                for child in parent.children
                    child.is_active = false
                end
                #Inherit status
                parent.status = first_status
                #Disconnect children to make parent a leaf
                empty!(parent.children)

                merges_count += 1
            else
            end
        end
    end
    println("  -> Merged $merges_count parents.")
end

"""
    generate_fine_to_coarse_map(mesh::CoarseMeshBuilder, partition::Tuple{Int, Int})

Generates a mapping between the uniform fine grid (defined by `partition`) and the current
coarse CoarseMeshBuilder.
"""
function generate_fine_to_coarse_map(mesh::CoarseMeshBuilder, partition::Tuple{Int, Int})
    nx, ny = partition

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
