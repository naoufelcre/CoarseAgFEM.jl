module Balancing

using ..QuadDefs
using ..CoreOps

export balance!

function balance!(mesh::CoarseMeshBuilder)
    println("  Phase 2: Balancing (Ripple Strategy)...")

    # Queue stores cells that need their neighborhood checked
    # Initial load: all active leaves
    queue = [n for n in mesh.all_nodes if n.is_active && is_leaf(n)]

    # We use a Set to avoid duplicate entries in queue (optimization)
    in_queue = Set{Int}()
    for n in queue
        push!(in_queue, n.id)
    end

    splits_performed = 0
    max_splits = 500000 # Safety brake

    while !isempty(queue)
        if splits_performed > max_splits
            println("  [ERROR] Balancing runaway! Exceeded $max_splits splits.")
            break
        end

        cell = popfirst!(queue)
        delete!(in_queue, cell.id)

        # If cell became inactive (split by someone else while in queue), skip
        if !cell.is_active || !is_leaf(cell) continue end

        # Check all 4 neighbors to see if there's a bigger cell in its neighbor, that's why 1.2 is enough.

        h = cell.size / 2.0
        dirs = [
            [0.0, h*1.2],  # N
            [0.0, -h*1.2], # S
            [h*1.2, 0.0],  # E
            [-h*1.2, 0.0]  # W
        ]

        for d in dirs
            p_check = cell.center .+ d

            # Boundary check
            if p_check[1] < 0.0 || p_check[1] > 1.0 || p_check[2] < 0.0 || p_check[2] > 1.0
                continue
            end

            neighbor = get_leaf_at(mesh.roots, p_check)

            # Rule: Neighbor cannot be >= 2 levels COARSER than cell.

            if (cell.level - neighbor.level) >= 2
                # Split Neighbor!
                split!(neighbor, mesh)
                splits_performed += 1

                # Add neighbor's NEW children to queue
                # They need to check their new environment
                for child in neighbor.children
                    if !(child.id in in_queue)
                        push!(queue, child)
                        push!(in_queue, child.id)
                    end
                end

            end
        end
    end

    println("  -> Balanced. Performed $splits_performed refinement ops.")
end

end # module
