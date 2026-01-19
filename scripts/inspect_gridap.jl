using Gridap
using GridapEmbedded

println("=== Gridap Inspection ===")
println("Path: ", pathof(Gridap))

# Check for keywords in names
keywords = ["hanging", "node", "conform", "junction", "constrain", "MultiPoint"]

function check_module(mod, name)
    println("\n--- $name ---")
    all_names = names(mod; all=true)
    
    found = false
    for n in all_names
        s = string(n)
        for k in keywords
            if occursin(k, lowercase(s))
                println("  Found symbol: $s")
                found = true
                break
            end
        end
    end
    if !found
        println("  No keywords found.")
    end
end

check_module(Gridap, "Gridap")
check_module(Gridap.Geometry, "Gridap.Geometry")
check_module(Gridap.FESpaces, "Gridap.FESpaces")
check_module(GridapEmbedded, "GridapEmbedded")

println("\n=== Testing H1 space on non-conforming mesh ===")
# Create a simple non-conforming 1D mesh? 2D hard to script without tools.
# Just checking documentation strings if possible?

println("Done.")
