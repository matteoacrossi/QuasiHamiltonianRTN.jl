# QuasiHamiltonianRTN

## Installation

    Pkg.clone("git://github.com/matteoacrossi/QuasiHamiltonianRTN.jl.git")

## Example

This example evaluates the 1D lattice

    using QuasiHamiltonianRTN
    using QuasiHamiltonianRTN.Lattice1D
    using PyPlot

    n = 10
    r0 = bloch_vector(localized_state(n,floor((n+1)/2)))

    t = linspace(0, 10, 100)

    # Noiseless qw
    H = Lattice1D.Hamiltonian(n, 0.)
    r = evolution(H, r0, collect(t), n; γ=0);

    imshow([real(diag(density_operator(r[:,i]))) for i in 1:length(t)], cmap="Blues")

    # Noisy qw with γ=1
    H = Lattice1D.Hamiltonian(n, 0.5)
    rn = evolution(H, r0, collect(t), n; γ=1.)

    imshow([real(diag(density_operator(rn[:,i]))) for i in 1:length(t)], cmap="Blues")


## TODO
* The MATLAB function expmv can evaluate for an equally spaced interval of t at
  once (using some optimization). Explore the possibility to do the same with `Expokit.jl`
