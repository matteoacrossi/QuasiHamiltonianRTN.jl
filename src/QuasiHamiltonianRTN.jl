module QuasiHamiltonianRTN
using Expokit # expmv!
include("SpecialUnitary.jl")
include("Utilities.jl")

include("Lattice1D.jl")
include("StarGraph.jl")

export quasiHamiltonian, quasiHamiltonian2, evolution, bloch_vector, density_operator, localized_state

"""
    Vnoise(n, γ)

Returns the generator ``V`` for the transition probabilities
between the various noise configurations for `n` fluctuators,
with switching rate ``\\gamma``
"""
function Vnoise(n, γ)
    V = sparse([-γ γ; γ -γ])

    function Vi(i)
        foldl(kron, [j == i ? V : speye(2) for j = 1:n])
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
#TODO: store a cache of structure constants (they don't use lots of space and
#      would save computational time)
function quasiHamiltonian(Hamiltonian::Function, RTN_number, γ)
    n = size(Hamiltonian(1), 1)
    Nc = 2^RTN_number
    f = structure_constants(n)
    λ = sun_generators(n)
    Hq = kron(-Vnoise(RTN_number,γ), speye(n^2-1))

    for i = 1 : Nc
        Hq += kron(D(RTN_number,i), igen(f, λ, Hamiltonian(i)))
    end
    return Hq
end

function quasiHamiltonian_old(Hamiltonian::Function, RTN_number, γ)
    n = size(Hamiltonian(1), 1)
    Nc = 2^RTN_number
    f = structure_constants(n)
    Hq = kron(-Vnoise(RTN_number,γ), speye(n^2-1))
    for i = 1 : Nc
        Hq += kron(D(n,i), igen_old(f, Hamiltonian(i)))
    end
    return Hq
end

# TODO: write more general version where we pass H0, Hnoise and the number of
# noise configurations
# function quasiHamiltonian(n, γ, ν)
#     f = structure_constants(n)
#     Hq = kron(-Vnoise(n,γ), speye(n^2-1))
#     Hq .+= kron(speye(2^n), igen(f,H0(n)))
#     Hq .+= ν * sum([kron(D(n,i), igen(f,Hn(n,i))) for i = 1:2^n])
#     return Hq
# end
#
# function quasiHamiltonian_old(n, γ, ν)
#     return kron(-Vnoise(n,γ), speye(n^2-1)) +
#         kron(speye(2^n), igen(H0(n))) +
#         ν * sum([kron(D(n,i), igen(Hnoise(n,i))) for i = 1:2^n])
# end

function D(n,i)
    D = spzeros(2^n,2^n)
    D[i,i] = 1
    return D
end

"""
    evolution(H, r0, t; [γ = 1])

Evaluates the dynamics of the system in presence of RTN noise with switching
rate γ, starting from an initial Bloch vector ``\\mathbf{r}_0`` of length
``N_Q = N^2-1``, where ``N`` is the dimension of the Hilbert space.

If t is an array of length ``N_T``, returns an ``N_T`` array of Bloch vectors
at each time instant.

## NOTE:
Due to the high cost of the construction of the quasi Hamiltonian, it is recommended
to evaluate simultaneously for all the time instants.
"""
#= TODO: Write a function where we pass the quasi Hamiltonian, so that there is no
        need to use the time vector. This can be useful e.g. in optimization problems
=#
function evolution(Hamiltonian::Function, r0, t, RTN_number; γ=1)
    n = size(Hamiltonian(1),1)

    Nq = length(r0)

    if n^2-1 != Nq
        throw(ArgumentError("Hamiltonian and r0 have incompatible dimensions"))
    end

    Hq = quasiHamiltonian(Hamiltonian, RTN_number, γ)

    Nc = 2^RTN_number

    rtnstate = ones(Nc)/sqrt(Nc)

    v0 = kron(rtnstate, r0)
    res = [
    begin
        tmp = zeros(Nq)
        v =  expmv(-ti, Hq, v0)
        #println(v)
        for k = 1:Nq
            for i = 1:Nc
                tmp[k] += rtnstate[i] * v[Nq*(i-1)+k]
            end
        end
        tmp
    end
    for ti in t]
    return res

#     y = kron(rtnstate, eye(Nq))
#     return y' * v
end

end # module
