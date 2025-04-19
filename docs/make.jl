using SP_IRK
using Documenter

DocMeta.setdocmeta!(SP_IRK, :DocTestSetup, :(using SP_IRK); recursive=true)

makedocs(;
    modules=[SP_IRK],
    authors="Fabio Durastante <fabio.durastante@unipi.it>, Mariarosa Mazza <mariarosa.mazza@uniroma2.it>",
    sitename="SP_IRK.jl",
    format=Documenter.HTML(;
        canonical="https://cirdans-home.github.io/SP_IRK.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cirdans-home/SP_IRK.jl",
    devbranch="main",
)
