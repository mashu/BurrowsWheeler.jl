using Documenter, BurrowsWheeler

makedocs(
    modules = [BurrowsWheeler],
    format = :html,
    checkdocs = :exports,
    sitename = "BurrowsWheeler.jl",
    pages = Any["index.md"],
    repo = "https://gitlab.com/matkad/burrowswheeler.jl/blob/{commit}{path}#{line}"
)

