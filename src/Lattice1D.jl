module Lattice1D
    export H0, Hnoise, Hamiltonian
    """
        H0(N; ϵ=0., ν=1.)

    Constructs the noiseless Hamiltonian for the 1D lattice with `N` nodes,
    self-energy ``\\epsilon`` and tunneling constant ``\\nu``.
    """
    function H0(n; ϵ=0., ν=1.)
        H = spdiagm((ϵ*ones(n), ν*ones(n-1), ν*ones(n-1)),(0,1,-1))
        if n > 2
            H[1,end] = H[end,1] = ν
        end
        return H
    end

    """
        Hnoise(N, noise_id)

    Returns the noisy part of the Hamiltonian for the 1D lattice with `N` nodes
    for a specific state of the noise fluctuators specified by `noise_id`.
    """
    function Hnoise(n, noise_id)
        if noise_id < 1 || noise_id > 2^n
            throw(ArgumentError("noise_id must be between 1 and 2^n"))
        end
        noise = 2 * digits(noise_id - 1, 2, n)-1
        H = spdiagm((noise[2:end], noise[2:end]), (-1,1))
        H[1, end] = H[end, 1] = noise[1]
        return H
    end

    """
        Hamiltonian(n, ν, ν0=1.)

    Returns an array of Hamiltonians of size `n` for all possible configurations
    of noise.

    # Arguments
        * `n::Integer`: the dimension of the Hilbert space
        * `ν::Float`: the strength of the noise
        * `ν0::Float`: The rate of the noiseless Hamiltonian (ν0 = 1.)
    """
    function Hamiltonian(n::Integer, ν::Float64, ν0=1.)
        return [H0(n;ν=ν0) + ν*Hnoise(n, id) for id = 1:2^n]
    end


end
