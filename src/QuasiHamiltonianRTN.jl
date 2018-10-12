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
with switching rate ``\\gamma``
"""
function Vnoise(n, γ)
    V = sparse([-γ γ; γ -γ])

    function Vi(i)
        foldl(kron, [j == i ? V : sparse(1.0I, 2, 2) for j = 1:n])
    end

    sum([Vi(i) for i = 1:n])
end

"""
    quasiHamiltonian(Hamiltonian, RTN_number, γ)

Returns the quasi-Hamiltonian ``H_q`` for the Hamiltonian specified by the
function `H` affected by a number `RTN_number` of RTN fluctuators with switching
rate γ.

The quasi-Hamiltonian is defined as in Eq. (20) of Joynt et al.

## Arguments
`Hamiltonian::Function`: H is a function of the type `Hamiltonian(i)`, where `i`
is an integer that goes from `1` to `2^RTN_number`, and returns a particular
configuration of the fluctuators.
"""
function quasiHamiltonian(Hamiltonian, γ)::SparseMatrixCSC{Float64,Int64}
    n = size(Hamiltonian[1], 1)
    Nc = length(Hamiltonian)
    RTN_number = Int(log2(Nc))
    f::Array{SparseArrays.SparseMatrixCSC{Float64,Int64},3} = structure_constants(n)
    λ::Array{SparseArrays.SparseMatrixCSC{Complex{Float64},Int64},1} = sun_generators(n)

    res::Array{SparseArrays.SparseMatrixCSC{Float64,Int64},1} = map(h -> igen(f,λ,h), Hamiltonian)
    Hq = blockdiag(res...) + kron(-Vnoise(RTN_number,γ), sparse(1.0I, n^2-1, n^2-1))

    return Hq
end

function D(n,i)
    D = spzeros(2^n,2^n)
    D[i,i] = 1
    return D
end

"""
    evolution(Hq, r0, t; [γ = 1])

Evaluates the dynamics of the system in presence of RTN noise with switching
rate γ, starting from an initial Bloch vector ``\\mathbf{r}_0`` of length
``N_Q = N^2-1``, where ``N`` is the dimension of the Hilbert space.

If t is an array of length ``N_T``, returns an ``N_T`` array of Bloch vectors
at each time instant.
"""
function evolution(Hq::SparseMatrixCSC, r0::Vector, t; γ=1)
    n = size(Hq, 1)
    # n = Nc * Nq
    Nq = length(r0)
    Nc = Int(n / Nq)

    rtnstate = ones(Nc)/sqrt(Nc)

    v0 = kron(rtnstate, r0)
    y = kron(rtnstate, sparse(1.0I, Nq, Nq))'
    res = y * ExpmV.expmv(-t, Hq, v0)

    return [res[:,i] for i = 1:size(res,2)]
end

function evolution(Hamiltonian::Array{SparseMatrixCSC{Float64,Int64},1},
            r0::Vector, t; γ=1)
    n = size(Hamiltonian[1], 1)
    Nc = length(Hamiltonian)
    Nq = length(r0)

    if n^2-1 != Nq
        throw(ArgumentError("Hamiltonian and r0 have incompatible dimensions"))
    end

    Hq = quasiHamiltonian(Hamiltonian, γ)

    return evolution(Hq, r0, t; γ=γ)
end


function evolution2(Hamiltonian, r0::Vector, t; γ=1)
    n = size(Hamiltonian[1], 1)
    Nc = length(Hamiltonian)
    Nq = length(r0)

    if n^2-1 != Nq
        throw(ArgumentError("Hamiltonian and r0 have incompatible dimensions"))
    end

    Hq = quasiHamiltonian(Hamiltonian, γ)

    rtnstate = ones(Nc)/sqrt(Nc)

    v0 = kron(rtnstate, r0)
    y = kron(rtnstate, sparse(1.0I, Nq, Nq))'
    res = [
    begin
        y * Expokit.expmv(-ti, Hq, v0)
    end
    for ti in t]
    return res
end

function evolution3(Hamiltonian, r0::Vector, t; γ=1)
    n = size(Hamiltonian[1], 1)
    Nc = length(Hamiltonian)
    Nq = length(r0)

    if n^2-1 != Nq
        throw(ArgumentError("Hamiltonian and r0 have incompatible dimensions"))
    end

    Hq = quasiHamiltonian(Hamiltonian, γ)

    rtnstate = ones(Nc)/sqrt(Nc)

    v0 = kron(rtnstate, r0)
    y = kron(rtnstate, sparse(1.0I, Nq, Nq))'
    res = [
    begin
        y * ExpmV.expmv(-ti, Hq, v0)
    end
    for ti in t]
    return hcat(res...)
end
end # module
