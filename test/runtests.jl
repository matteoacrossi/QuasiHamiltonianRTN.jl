using QuasiHamiltonianRTN
using QuasiHamiltonianRTN.Lattice1D
using Combinatorics
using SparseArrays
using LinearAlgebra

include("../src/SpecialUnitary.jl")
include("../src/Utilities.jl")

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@testset "SU(N) generators" begin
    @testset "N = $N" for N = 20:21
        @testset "trless" begin
            @test mapreduce(A -> isapprox(tr(A), 0., atol=eps(Float64) * N), &, sun_generators(N))
        end
        # Hermitian
        @testset "Hermitian" begin
            @test mapreduce(A -> A ≈ A', &, sun_generators(N))
        end
        # Orthogonal
        @testset "Orthogonal" begin
            @test let λ = sun_generators(N)
                [tr(λ[i] * λ[j]) for i = 1:length(λ), j = 1:length(λ)] ≈ 2*I
             end
        end
    end
end

@testset "SU(N) structure constants" begin
    @testset "SU(2) gives the Levi-Civita tensor" begin
        @test reduce(&, [structure_constants(2)[i][j,k] == levicivita([i,j,k]) for i=1:3,j=1:3,k=1:3])
    end

    @testset "[λ[i], λ[j]] = i ∑ f[i][j,k] λ[k] for SU($n)" for n = 7:9
        let λ = sun_generators(n), f = structure_constants(n)
            @test reduce(&,[commutator(λ[i], λ[j]) ≈ 2*im*sum([f[i][j,k]*λ[k] for k = 1:n^2-1])  for i = 1 : n^2-1, j = 1:n^2-1])
        end
    end
end

@testset "Bloch vector" begin
    @testset "Back and forth" begin
        for n = 2:10
            ρ = random_density_matrix(n)
            @test Matrix(density_operator(bloch_vector(ρ))) ≈ ρ
        end
    end
end

@testset "Unitary evolution check" begin
    let n = 4
        ρ = random_density_matrix(n)
        ρ1 = unitary_op(Lattice1D.H0(n), 1.)*ρ*unitary_op(Lattice1D.H0(n), 1.)'
        ρ2 = Matrix(density_operator(exp(- Matrix(igen(Lattice1D.H0(n)))) * bloch_vector(ρ)))
        @test ρ1 ≈ ρ2
    end
end

@testset "Noiseless evolution check" begin
    let tv = range(0, stop=5, length=10), n = 5, ρ0 = localized_state(n, 3)
        rtU = [bloch_vector(unitary_op(Lattice1D.H0(n), t) * ρ0 * unitary_op(Lattice1D.H0(n), t)') for t in tv]
        rtE = evolution(Lattice1D.Hamiltonian(n, .0), bloch_vector(ρ0), tv, γ=0.)
        @test rtU ≈ rtE

        Hq = quasiHamiltonian(Lattice1D.Hamiltonian(n, .0,), 0.)
        rtEHq = evolution(Hq, bloch_vector(ρ0), tv, γ=0.)

        @test rtEHq ≈ rtE
    end
end
