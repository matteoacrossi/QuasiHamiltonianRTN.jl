using Documenter
using QuasiHamiltonianRTN

makedocs(
    sitename = "QuasiHamiltonianRTN",
    format = :html,
    modules = [QuasiHamiltonianRTN]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/matteoacrossi/QuasiHamiltonianRTN.jl.git",
)
