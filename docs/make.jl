using BurrowsWheeler
using Documenter

makedocs(;
    modules=[BurrowsWheeler],
    authors="Mateusz Kaduk <mateusz.kaduk@ki.se> and contributors",
    sitename="BurrowsWheeler.jl",
    format=Documenter.HTML(;
        canonical="https://mashu.github.io/BurrowsWheeler.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Overview" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mashu/BurrowsWheeler.jl",
    devbranch="main",
    versions=nothing,
)

