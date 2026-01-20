module Paving

using ..QuadDefs
using ..CoreOps

export pave_mesh

function analyze_faces(leaf::QuadNode, roots::Vector{QuadNode})
    faces = Dict(:N=>false, :E=>false, :S=>false, :W=>false)
    hanging_nodes = []

    h = leaf.size / 2.0
    c = leaf.center

    dirs = [
        (:N, [0.0, h*1.1], [0.0, h]),
        (:S, [0.0, -h*1.1], [0.0, -h]),
        (:E, [h*1.1, 0.0], [h, 0.0]),
        (:W, [-h*1.1, 0.0], [-h, 0.0])
    ]
    face_map = Dict{Symbol, Vector{Float64}}()

    for (name, probe_offset, node_offset) in dirs
        p = c .+ probe_offset
        if p[1] < 0 || p[1] > 1 || p[2] < 0 || p[2] > 1 continue end

        neigh = get_leaf_at(roots, p)
        if neigh.level > leaf.level
            hn = c .+ node_offset
            push!(hanging_nodes, hn)
            face_map[name] = hn
        end
    end
    return face_map, length(hanging_nodes)
end

function get_face_mask(leaf::QuadNode, roots::Vector{QuadNode})
    mask = 0

    h = leaf.size / 2.0
    c = leaf.center

    # Epsilon for robust probing
    eps = h * 0.01

    #BitValue allows the code to compress the state of all 4 neighbors into a single unique integer ID

    probes = [
        # North (Bit 1): P1=Just Left of N-mid, P2=Just Right of N-mid
        (:N, 1, [-eps, h+eps], [eps, h+eps]),

        # South (Bit 2): P1=Just Left of S-mid, P2=Just Right of S-mid
        (:S, 2, [-eps, -h-eps], [eps, -h-eps]),

        # East (Bit 4): P1=Just Below E-mid, P2=Just Above E-mid
        (:E, 4, [h+eps, -eps], [h+eps, eps]),

        # West (Bit 8): P1=Just Below W-mid, P2=Just Above W-mid
        (:W, 8, [-h-eps, -eps], [-h-eps, eps])
    ]

    # Hanging nodes locations for construction
    hn_dict = Dict{Symbol, Vector{Float64}}()

    for (name, bit, off1, off2) in probes
        p1 = c .+ off1
        p2 = c .+ off2

        # Bounds check
        in_bounds1 = (p1[1] >= 0 && p1[1] <= 1 && p1[2] >= 0 && p1[2] <= 1)
        in_bounds2 = (p2[1] >= 0 && p2[1] <= 1 && p2[2] >= 0 && p2[2] <= 1)

        if !in_bounds1 || !in_bounds2 continue end

        n1 = get_leaf_at(roots, p1)
        n2 = get_leaf_at(roots, p2)

        # If n1 != n2, it means the neighbor across this face is SPLIT.
        # This implies a hanging node exists at the face midpoint.
        if n1.id != n2.id
            mask |= bit

            # DETERMINISTIC SNAPPING
            # We use the neighbor 'n1' found by the first probe to determine the exact coordinate of the hanging node.

            h1 = n1.size / 2.0
            c1 = n1.center

            hn_pos = [0.0, 0.0]

            if name == :N
                # P1 checked NW direction (y > center). n1 is the "Left" neighbor on top.
                # Its SE corner touches the hanging node.
                hn_pos = c1 .+ [h1, -h1]
                hn_dict[:N] = hn_pos

            elseif name == :S
                # P1 checked SW direction (y < center). n1 is the "Left" neighbor on bottom.
                # Its NE corner touches the hanging node.
                hn_pos = c1 .+ [h1, h1]
                hn_dict[:S] = hn_pos

            elseif name == :E
                # P1 checked SE direction (x > center, y < center). n1 is the "Bottom" neighbor on right.
                # Its NW corner touches the hanging node.
                hn_pos = c1 .+ [-h1, h1]
                hn_dict[:E] = hn_pos

            elseif name == :W
                # P1 checked SW direction (x < center, y < center). n1 is the "Bottom" neighbor on left.
                # Its NE corner touches the hanging node.
                hn_pos = c1 .+ [h1, h1]
                hn_dict[:W] = hn_pos
            end
        end
    end

    return mask, hn_dict
end

function propagate_splits!(mesh::CoarseMeshBuilder)
    println("  Phase 2b: Propagating Splits (Odd -> Even Hanging Nodes)...")

    # 1. Initialize Queue with all active leaves
    queue = Set{Int}()
    leaves = [n for n in mesh.all_nodes if n.is_active && is_leaf(n)]
    for n in leaves
        push!(queue, n.id)
    end

    q_vec = leaves
    q_set = Set{Int}(n.id for n in leaves)

    splits = 0
    max_splits = 1000000 # Increased for high-res meshes

    while !isempty(q_vec)
        if splits > max_splits
            println("  [WARNING] Split propagation limit reached ($max_splits)!")
            break
        end

        leaf = popfirst!(q_vec)
        delete!(q_set, leaf.id)

        if !leaf.is_active || !is_leaf(leaf) continue end

        mask, _ = get_face_mask(leaf, mesh.roots)
        n_hanging = count_ones(mask)

        # Rule: If Odd (1 or 3), SPLIT.
        if n_hanging == 1 || n_hanging == 3
            split!(leaf, mesh)
            splits += 1

            # Add new children to queue
            for child in leaf.children
                if parse(Int, string(child.id)) ∉ q_set # optimization? id is Int.
                    push!(q_vec, child)
                    push!(q_set, child.id)
                end
            end

            #We changed coordinates of the interface.
            # Neighbors now have potentially NEW hanging nodes (the children of this split).
            # We must re-evaluate neighbors.
            h = leaf.size / 2.0
            dirs = ([0, h*1.1], [0, -h*1.1], [h*1.1, 0], [-h*1.1, 0])
            for d in dirs
                nb = get_leaf_at(mesh.roots, leaf.center .+ d)
                if nb.id != leaf.id && nb.is_active && is_leaf(nb)
                    if !(nb.id in q_set)
                        push!(q_vec, nb)
                        push!(q_set, nb.id)
                    end
                end
            end
        end
    end
    println("  -> Propagated $splits additional splits.")
end

function pave_mesh(mesh::CoarseMeshBuilder)
    # Pre-Process: Eliminate Odd Hanging Nodes
    # propagate_splits!(mesh)

    println("  Phase 3: Template Paving (All-Quad Patterns)...")
    elements = QuadElement[]
    leaves = [n for n in mesh.all_nodes if n.is_active && is_leaf(n)]

    for leaf in leaves
        mask, hn = get_face_mask(leaf, mesh.roots)

        h = leaf.size / 2.0
        c = leaf.center

        # Corners
        nw = c .+ [-h, h]
        ne = c .+ [ h, h]
        sw = c .+ [-h,-h]
        se = c .+ [ h,-h]

        col = "cornflowerblue"
        col_trans = "orange"

        # --- Templates ---

        if mask == 0
            # Clean Quad
            push!(elements, QuadElement(leaf.id, [sw, se, ne, nw], col))

        # --- Tunnel Cases (2 Hanging Nodes, Opposite) ---
        elseif mask == 3 # N(1) + S(2)
            # Vertical Split
            # Left Quad: SW, S_mid, N_mid, NW
            q1 = [sw, hn[:S], hn[:N], nw]
            # Right Quad: S_mid, SE, NE, N_mid
            q2 = [hn[:S], se, ne, hn[:N]]
            push!(elements, QuadElement(leaf.id, q1, col_trans))
            push!(elements, QuadElement(leaf.id, q2, col_trans))

        elseif mask == 12 # E(4) + W(8)
            # Horizontal Split
            # Bottom Quad: SW, SE, E_mid, W_mid
            q1 = [sw, se, hn[:E], hn[:W]]
            # Top Quad: W_mid, E_mid, NE, NW
            q2 = [hn[:W], hn[:E], ne, nw]
            push!(elements, QuadElement(leaf.id, q1, col_trans))
            push!(elements, QuadElement(leaf.id, q2, col_trans))

        # --- Corner / Kite Cases (2 Hanging Nodes, Adjacent) ---
        # The logic: 3 Quads.
        # We need an internal node 'm'. Paving.jl usually doesn't create internal nodes?
        # Actually, we can just compute leaf.center as 'm'.

        # S(2) + E(4) = 6
        elseif mask == 6

            # Mask 6 (S+E)
            q1 = [sw, hn[:S], c, nw]
            q2 = [hn[:S], se, hn[:E], c]
            q3 = [c, hn[:E], ne, nw] # Note winding: c->e->ne->nw
            push!(elements, QuadElement(leaf.id, q1, col_trans))
            push!(elements, QuadElement(leaf.id, q2, col_trans))
            push!(elements, QuadElement(leaf.id, q3, col_trans))

        # S(2) + W(8) = 10
        elseif mask == 10

            q1 = [hn[:W], sw, hn[:S], c]
            q2 = [hn[:S], se, ne, c]
            q3 = [c, ne, nw, hn[:W]]
            push!(elements, QuadElement(leaf.id, q1, col_trans))
            push!(elements, QuadElement(leaf.id, q2, col_trans))
            push!(elements, QuadElement(leaf.id, q3, col_trans))

        # N(1) + E(4) = 5
        elseif mask == 5

            q1 = [c, hn[:E], ne, hn[:N]]
            q2 = [sw, se, hn[:E], c]
            q3 = [sw, c, hn[:N], nw]
            push!(elements, QuadElement(leaf.id, q1, col_trans))
            push!(elements, QuadElement(leaf.id, q2, col_trans))
            push!(elements, QuadElement(leaf.id, q3, col_trans))

        # N(1) + W(8) = 9
        elseif mask == 9

            q1 = [hn[:W], c, hn[:N], nw] # NW Corner
            q2 = [c, se, ne, hn[:N]]     # East Block
            q3 = [sw, se, c, hn[:W]]     # South Block
            push!(elements, QuadElement(leaf.id, q1, col_trans))
            push!(elements, QuadElement(leaf.id, q2, col_trans))
            push!(elements, QuadElement(leaf.id, q3, col_trans))

        elseif mask == 15
            # 4 Hanging Nodes -> Use Cross Split (4 Quads)
            q1 = [sw, hn[:S], c, hn[:W]]
            q2 = [hn[:S], se, hn[:E], c]
            q3 = [c, hn[:E], ne, hn[:N]]
            q4 = [hn[:W], c, hn[:N], nw]
            push!(elements, QuadElement(leaf.id, q1, col_trans))
            push!(elements, QuadElement(leaf.id, q2, col_trans))
            push!(elements, QuadElement(leaf.id, q3, col_trans))
            push!(elements, QuadElement(leaf.id, q4, col_trans))

        else
            # ODD MASKS (1, 2, 4, 8, 7, 11, 13, 14) -> SAFE TRIANGLE FAN

            ring = []

            # Start SW (Bottom-Left)
            push!(ring, sw)

            # South Face
            if (mask & 2) != 0 push!(ring, hn[:S]) end

            # SE Corner
            push!(ring, se)

            # East Face
            if (mask & 4) != 0 push!(ring, hn[:E]) end

            # NE Corner
            push!(ring, ne)

            # North Face
            if (mask & 1) != 0 push!(ring, hn[:N]) end

            # NW Corner
            push!(ring, nw)

            # West Face
            if (mask & 8) != 0 push!(ring, hn[:W]) end

            # Generate Triangles
            col_fan = "yellow"
            for i in 1:length(ring)
                p1 = ring[i]
                p2 = ring[mod1(i+1, length(ring))]

                # Connectivity: c -> p1 -> p2 (CCW)
                # Valid because Center is inside, and p1->p2 is CCW on boundary.
                push!(elements, QuadElement(leaf.id, [c, p1, p2], col_fan))
            end
        end
    end

    return elements
end



end # module
