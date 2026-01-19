using CoarseAgFEM
using Documenter

makedocs(;
    modules=[CoarseAgFEM],
    authors="Naoufel Cresson",
    repo="https://github.com/Naoufel/CoarseAgFEM.jl/blob/{commit}{path}#{line}",
    sitename="CoarseAgFEM.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Naoufel.github.io/CoarseAgFEM.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Naoufel/CoarseAgFEM.jl",
    devbranch="main",
)
