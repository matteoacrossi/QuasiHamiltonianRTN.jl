module StarGraph
    export H0, Hnoise, Hamiltonian
    """
        H0(N; ν=1.)

    Constructs the noiseless Hamiltonian for the star graph with ``N`` nodes
    """
    function H0(N; ν=1.)
        H = spzeros(N,N)
        H[1,2:end] = ν
        H = H + H'
        H = spdiagm(0 => sum(H,1)[:]) - H
        return H
    end

    """
        Hnoise(N, noise_id)

    Returns the noisy part of the Hamiltonian for the 1D lattice with `N` nodes
    for a specific state of the noise fluctuators specified by `noise_id`.
    """
    function Hnoise(N, noise_id)
        if noise_id < 0 || noise_id > 2^(N-1)
            throw(ArgumentError("noise_id must be between 1 and 2^(N-1)"))
        end

        if noise_id ==0
            return spzeros(N,N)
        end

        noise = 2 * digits(noise_id - 1, base=2, pad=N-1) - 1

        H = spzeros(N,N)
        H[1, 2:end] = noise
        H = H + H'
        H = spdiagm(0 => sum(H,1)[:]) - H
        return H
    end

    """
        Hamiltonian(N, ν, ν0=1.)

    Returns an array of Hamiltonians for all possible configurations of the noise

    # Arguments
        * `N::Integer`: the dimension of the Hilbert space
        * `ν::Float`: the strength of the noise
        * `ν0::Float`: The rate of the noiseless Hamiltonian (ν0 = 1.)
    """
    function Hamiltonian(N::Integer, ν::Float64, ν0=1.)
        return [H0(N;ν=ν0) + ν*Hnoise(N, id) for id = 1:2^(N-1)]
    end

    function H0_spatial_search(N; ν=1., target::Int=1)
        H = spzeros(N,N)
        H[1,2:end] = ν
        H = H + H'
        H = spdiagm(0 => sum(H,1)[:]) - H
        H[target,target] = -1.
        return H
    end

    function Hamiltonian_spatial_search(N::Integer, ν::Float64; ν0=.3, target::Int=1)
        return [H0_spatial_search(N;ν=ν0) + ν*Hnoise(N, id) for id = 1:2^(N-1)]
    end
end
