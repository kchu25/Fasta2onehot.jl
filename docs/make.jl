using Fasta2onehot
using Documenter

DocMeta.setdocmeta!(Fasta2onehot, :DocTestSetup, :(using Fasta2onehot); recursive=true)

makedocs(;
    modules=[Fasta2onehot],
    authors="Shane Kuei-Hsien Chu (skchu@wustl.edu)",
    sitename="Fasta2onehot.jl",
    format=Documenter.HTML(;
        canonical="https://kchu25.github.io/Fasta2onehot.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kchu25/Fasta2onehot.jl",
    devbranch="main",
)
