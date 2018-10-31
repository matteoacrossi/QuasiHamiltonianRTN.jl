module QuasiHamiltonianRTN
using Expokit # expmv
using ExpmV
using SparseArrays
using LinearAlgebra

include("SpecialUnitary.jl")
include("Utilities.jl")

include("Lattice1D.jl")
include("StarGraph.jl")

export quasiHamiltonian, evolution, bloch_vector, density_operator, localized_state

"""
    Vnoise(n, γ)

Returns the generator ``V`` for the transition probabilities
between the various noise configurations for `n` fluctuators,
with switching rate ``\\gamma``.

It is defined as

``\\mathbf{V} = \\sum_{i=1}^N \\mathbf{V}_\\gamma^{(i)}``

with ``\\mathbf{V}_\\gamma^{(i)} = \\mathbf{I}_2^{\\otimes i-1} \\otimes \\mathbf{V}_\\gamma \\otimes \\mathbf{I}_2^{\\otimes N-i}``.
"""
function Vnoise(n, γ)
    V = sparse([-γ γ; γ -γ])

    function Vi(i)
        foldl(kron, [j == i ? V : sparse(1.0I, 2, 2) for j = 1:n])
    end

    sum([Vi(i) for i = 1:n])
end

"""
    quasiHamiltonian(Hamiltonian, γ)

Returns the quasi-Hamiltonian ``H_q`` for the Hamiltonian specified by the
array of matrices `H`, representing all the possible configurations of the
Hamiltonian, with switching rate γ.

The quasi-Hamiltonian is defined as in Eq. (20) of Joynt et al.
This function returns ``-i H_q``, in order to deal with real matrices instead of
imaginary ones.

## Arguments

* `Hamiltonian` is an array of matrices, each of them representing the Hamiltonian
with a different configuration of the noise.
* `γ` is the switching rate of the noise.
"""
function quasiHamiltonian(Hamiltonian, γ)::SparseMatrixCSC{Float64,Int64}
    n = size(Hamiltonian[1], 1)
    Nc = length(Hamiltonian)
    RTN_number = Int(log2(Nc))
    λ::Array{SparseArrays.SparseMatrixCSC{Complex{Float64},Int64},1} = sun_generators(n)
    f::Array{SparseArrays.SparseMatrixCSC{Float64,Int64},3} = structure_constants(λ)

    res::Array{SparseArrays.SparseMatrixCSC{Float64,Int64},1} = map(h -> igen(f,λ,h), Hamiltonian)
    Hq = blockdiag(res...) + kron(-Vnoise(RTN_number,γ), sparse(1.0I, n^2-1, n^2-1))

    return Hq
end

"""
    evolution(Hq, r0, t; [expmpkg=:ExpmV])

Evaluates the state of the system with quasi-Hamiltonian matrix `Hq`, starting
from an initial Bloch vector `r0`, at time `t`.

If `t` is a `StepRangeLen` object of length ```N_t``, return a ``N_r \\times N_t``
array containing the Bloch vector of the system at each time step as columns.
"""
function evolution(Hq::SparseMatrixCSC, r0::Vector, t; expmpkg=:ExpmV)
    @assert expmpkg ∈ (:ExpmV, :Expokit) "expmpkg must be either ExpmV or Expokit"

    n = size(Hq, 1)
    # n = Nc * Nq
    Nq = length(r0)
    Nc = Int(n / Nq)

    rtnstate = ones(Nc)/sqrt(Nc)

    v0 = kron(rtnstate, r0)
    y = kron(rtnstate, sparse(1.0I, Nq, Nq))'

    if expmpkg == :ExpmV
        res = ExpmV.expmv(-t, Hq, v0)

    elseif expmpkg == :Expokit
        res = Array{Float64}(undef, n, length(t))
        for (i, ti) in enumerate(t)
            res[:, i] = Expokit.expmv(-ti, Hq, v0)
        end

    end

    return y * res
end

function evolution(Hamiltonian::Array{SparseMatrixCSC{Float64,Int64},1},
            r0::Vector, t; γ=1, expmpkg=:ExpmV)
    n = size(Hamiltonian[1], 1)
    Nc = length(Hamiltonian)
    Nq = length(r0)

    if n^2-1 != Nq
        throw(ArgumentError("Hamiltonian and r0 have incompatible dimensions"))
    end

    Hq = quasiHamiltonian(Hamiltonian, γ)

    return evolution(Hq, r0, t; expmpkg=expmpkg)
end

end # module
