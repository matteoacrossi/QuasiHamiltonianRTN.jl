# QuasiHamiltonianRTN
[![Build Status](https://travis-ci.org/matteoacrossi/QuasiHamiltonianRTN.jl.svg?branch=master)](https://travis-ci.org/matteoacrossi/QuasiHamiltonianRTN.jl) [![codecov](https://codecov.io/gh/matteoacrossi/QuasiHamiltonianRTN.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/matteoacrossi/QuasiHamiltonianRTN.jl)

## Installation

    Pkg.clone("git://github.com/matteoacrossi/QuasiHamiltonianRTN.jl.git")

## Example

This example evaluates the 1D lattice

```julia
using QuasiHamiltonianRTN
using QuasiHamiltonianRTN.Lattice1D
using PyPlot

n = 10
r0 = bloch_vector(localized_state(n,floor((n+1)/2)))

t = range(0, stop=10, length=100)

# Noiseless qw
H = Lattice1D.Hamiltonian(n, 0.)
r = evolution(H, r0, collect(t), n; γ=0);

imshow([real(diag(density_operator(r[:,i]))) for i in 1:length(t)], cmap="Blues")

# Noisy qw with γ=1
H = Lattice1D.Hamiltonian(n, 0.5)
rn = evolution(H, r0, collect(t), n; γ=1.)

imshow([real(diag(density_operator(rn[:,i]))) for i in 1:length(t)], cmap="Blues")
```

## TODO
* The MATLAB function [`expmv`](https://github.com/higham/expmv) can evaluate ``exp(t A)`` for an equally spaced interval of t at
  once (using some optimization). Explore the possibility to do the same with `Expokit.jl`, or, most probably, extend the code in [ExpmV.jl](https://github.com/marcusps/ExpmV.jl) which is the direct port of `expmv`.
