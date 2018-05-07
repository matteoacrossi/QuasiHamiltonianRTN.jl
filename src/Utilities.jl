"""
	random_density(n)

Returns a random ``N × N`` density matrix, distributed according to the
Hilbert-Schmidt measure.

NOTE: No guarantee about the quality of the algorithm
For more details see [V Al Osipov et al 2010 J. Phys. A: Math. Theor. 43 055302]
(https://doi.org/10.1088/1751-8113/43/5/055302)
"""
function random_density_matrix(n)
    A = randn(n,n) + 1im * randn(n,n)
    rho = A * A'
    return rho / trace(rho)
end

"""
	unitary_op(H, t)

Returns the unitary evolution operator ``U = \\exp[-i H t]``.
"""
function unitary_op(H, t)
   return expm(-1im * full(H) * t)
end

"""
    commutator(A, B)

Returns the commutator [A,B] = A * B - B * A between two matrices
"""
function commutator(a, b)
	a * b - b * a
end

"""
	localized_state(N, i)

Returns the density operator for a state ``|i⟩`` in an ``N``-dimensional Hilbert
space
"""
function localized_state(n, i)
   sparse([i,n],[i,n],[1,0])
end

"""
    probability_current(ρ, A, c)

Returns the probability current for the node ``|c ⟩``,
defined as the expected value of the operator

``C_c = i\\sum_{⟨j,c⟩}(|j⟩⟨c| - |c⟩⟨j|)``,

where the sum is over the neighbors of the node ``|c⟩``.

```math
 \\text{Tr} [C_c ρ] = - 2 ∑_{j≠c}ℑ(A_{cj}ρ_{cj})
```
"""
function probability_current(ρ::AbstractArray{Complex128,2}, A::Matrix, c::Integer)
    return - 2 * imag(A[c,:] * ρ[c,:]) # NOTE: ρ[c,c] is real
end
