using SequentialMeasurements
using Documenter

DocMeta.setdocmeta!(SequentialMeasurements, :DocTestSetup, :(using SequentialMeasurements); recursive=true)

makedocs(;
    modules=[SequentialMeasurements],
    authors="Diego N. Bernal-Garcia",
    sitename="SequentialMeasurements.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
