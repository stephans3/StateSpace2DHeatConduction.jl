using StateSpace2DHeatConduction
using Documenter

DocMeta.setdocmeta!(StateSpace2DHeatConduction, :DocTestSetup, :(using StateSpace2DHeatConduction); recursive=true)

makedocs(;
    modules=[StateSpace2DHeatConduction],
    authors="Stephan Scholz",
    sitename="StateSpace2DHeatConduction.jl",
    format=Documenter.HTML(;
        canonical="https://stephans3.github.io/StateSpace2DHeatConduction.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/stephans3/StateSpace2DHeatConduction.jl",
    devbranch="main",
)
