"""
    sun_generators(n)

Returns an array of sparse ``n \\times n`` generators of ``SU(N)``
that are Hermitian and traceless.
"""
function sun_generators(n)
    matrices = [spzeros(Complex128, n, n) for i = 1 : n^2-1]
    i = 1
    for k = 2:n
        for j = 1:k-1
            matrices[i][j, k] = matrices[i][k, j] = 1.
            i += 1
            matrices[i][j, k] -= (matrices[i][k, j] = 1im)
            i += 1
        end
    end

    for l = 1:n-1
        for j = 1:l
            matrices[i][j, j] = sqrt(2/(l*(l+1)))
        end
        matrices[i][l+1, l+1] = -sqrt(2*l/(l+1))
        i += 1
    end
    return matrices
end

"""
    bloch_vector(ρ)

Returns the generalized Bloch vector ``\\mathbf{r}`` for the density operator ρ.

If ρ is a ``N × N`` density operator, ``\\mathbf{r}`` will be a ``N^2 - 1``
vector of real values. We follow the formula
``ρ = \\frac{1}{N} (I + \\sqrt{N} \\mathbf{r} \\cdot \\mathbf{λ})``
"""
function bloch_vector(ρ)
   n = size(ρ)[1]
   λ = sun_generators(n)
   sqrt(n)/2*real([trace(λ[i] * ρ) for i = 1: length(λ)])
end

"""
    density_operator(r)

Returns the density operator ρ corresponding to the Bloch vector `r`.

`vec` must be a vector of ``N^2-1`` real values. We follow the formula
``ρ = \\frac{1}{N} (I + \\sqrt{N} \\mathbf{r} \\cdot \\mathbf{λ})``
"""
function density_operator(vec::Array{Float64})
   n = length(vec)
   Nq = Int(floor(sqrt(n+1))) # Dimension of the corresponding Hilbert space
   λ = sun_generators(Nq) # Obtain the generators of the SU(N) Hilbert space

   1/Nq*(I + sqrt(Nq) * sum([vec[i] * λ[i] for i = 1:n]))
end

# We return im * gen so that the matrix is real (more efficient)


"""
    igen(H)

Returns ``i G(H)``, where G is defined in Eq. (11) of Joynt et al.
"""
function igen(H)
    n = size(H,1)
    λ = sun_generators(n)
    result = spzeros(n^2-1, n^2-1)
    for i = 1:n^2-1
        for j = i+1:n^2-1
            tmp = imag(0.5 * trace(commutator(λ[i],λ[j])*H))
            if tmp != 0
                @inbounds result[i,j] = tmp
                @inbounds result[j,i] = -tmp
            end
        end
    end
    return result
end

"""

    igen(f, H)

Returns ``i G(H)``, where G is defined in Eq. (11) of Joynt et al., with the
structure constants tensor `f` that can be evaluated with `structure_constants(n)`

When providing `f`, the calculations are much faster (especially if many
generators have to be evaluated).
"""
function igen(f, H)
    n = size(H,1)
    λ = sun_generators(n)
    result = spzeros(n^2-1, n^2-1)
    ak = real([trace(H*lambda) for lambda in λ])
    for i = 1:n^2-1
        @inbounds result[i,:]= f[i] * ak
    end
    return result
end

"""
    structure_constants(N)

Returns the totally antisymmetric tensor ``i f_{ijk}`` of the structure constants
of ``SU(N)``, defined as ``f_{ijk} = - 2 i \\text{Tr}([λ_i,λ_j]λ_k)``.

Notice that the true definition has a factor ``-i``, which we avoid to store only
a real value.
"""
function structure_constants(n)
    lambda = sun_generators(n)
    f = cat(3, [spzeros(n^2-1,n^2-1) for i = 1:n^2-1])
    cacca = 0.
    for i = 1:n^2-1
        for j = 1:n^2-1
            for k = j+1:n^2-1
                @inbounds tmp = 1/4*imag(trace(commutator(lambda[i],lambda[j])*lambda[k]))
                @inbounds f[i][j,k] = tmp
                @inbounds f[i][k,j] = -tmp
            end
        end
    end
    return f
end