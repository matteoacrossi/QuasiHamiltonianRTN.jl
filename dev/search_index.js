var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "QuasiHamiltonianRTN.jl",
    "title": "QuasiHamiltonianRTN.jl",
    "category": "page",
    "text": ""
},

{
    "location": "#QuasiHamiltonianRTN.jl-1",
    "page": "QuasiHamiltonianRTN.jl",
    "title": "QuasiHamiltonianRTN.jl",
    "category": "section",
    "text": "This package allows for the numerical solution of the dynamics of a quantum system affected by classical random telegraph noise (RTN), using the quasi-Hamiltonian method introduced in Joynt et al. (2011).NOTE: The algorithm complexity is linear in the number of noise states, meaning that it is exponential for independent fluctuators."
},

{
    "location": "#Model-1",
    "page": "QuasiHamiltonianRTN.jl",
    "title": "Model",
    "category": "section",
    "text": "Consider a quantum system described by the density operator rho(t), with an Hamiltonian that is function of one or more stochastic processes g_i(t). The dynamics of the system is described by the dynamical mapbeginequation\nlabeleqmap\nrho(t) = leftlangle mathcalT e^-i Hg_i(t) t  rho(0)  mathcalT e^i Hg_i(t) t rightrangle_g_i(t)\nendequationWhen the number of possible configurations of the classical environment is small enough, we can employ the exact method described in Details on the method.A typical case of use is RTN: a stochastic noise which can have two values, for example (pm 1), and jumps between them with a certain switching rate gamma.(Image: )The package QuasiHamiltonianRTN allows for obtaining a numerically exact solution of the dynamics without recurring to Montecarlo simulations of the equation \\eqref{eq:map}."
},

{
    "location": "#Usage-1",
    "page": "QuasiHamiltonianRTN.jl",
    "title": "Usage",
    "category": "section",
    "text": "To use the package, first construct an array of matrices, in which each element corresponds to the Hamiltonian with a particular realization of the noise sources.For instance, if there are two noise sources, the array will be composed of four Hamiltonians.Once the vector has been defined, say H, we generate the quasi-Hamiltonian matric"
},

{
    "location": "#QuasiHamiltonianRTN.evolution",
    "page": "QuasiHamiltonianRTN.jl",
    "title": "QuasiHamiltonianRTN.evolution",
    "category": "function",
    "text": "evolution(Hq, r0, t; [expmpkg=:ExpmV])\n\nEvaluates the state of the system with quasi-Hamiltonian matrix Hq, starting from an initial Bloch vector r0, at time t.\n\nIf t is a StepRangeLen object of length `N_t, return a N_r times N_t array containing the Bloch vector of the system at each time step as columns.\n\n\n\n\n\n"
},

{
    "location": "#QuasiHamiltonianRTN.igen",
    "page": "QuasiHamiltonianRTN.jl",
    "title": "QuasiHamiltonianRTN.igen",
    "category": "function",
    "text": "igen(H)\n\nReturns i G(H), where G is defined in Eq. (11) of Joynt et al.\n\n\n\n\n\nigen(f, λ, H)\n\nReturns i G(H), where G is defined in Eq. (11) of Joynt et al., with the structure constants tensor f that can be evaluated with structure_constants(λ), and the array of SU(N) generators λ (sun_generators(n)).\n\nWhen providing f, the calculations are much faster (especially if many generators have to be evaluated).\n\n\n\n\n\n"
},

{
    "location": "#QuasiHamiltonianRTN.bloch_vector",
    "page": "QuasiHamiltonianRTN.jl",
    "title": "QuasiHamiltonianRTN.bloch_vector",
    "category": "function",
    "text": "bloch_vector(ρ)\n\nReturns the generalized Bloch vector mathbfr for the density operator ρ.\n\nIf ρ is a N  N density operator, mathbfr will be a N^2 - 1 vector of real values. We follow the formula ρ = frac1N (I + sqrtN mathbfr cdot mathbfλ)\n\n\n\n\n\n"
},

{
    "location": "#QuasiHamiltonianRTN.quasiHamiltonian",
    "page": "QuasiHamiltonianRTN.jl",
    "title": "QuasiHamiltonianRTN.quasiHamiltonian",
    "category": "function",
    "text": "quasiHamiltonian(Hamiltonian, γ)\n\nReturns the quasi-Hamiltonian H_q for the Hamiltonian specified by the array of matrices H, representing all the possible configurations of the Hamiltonian, with switching rate γ.\n\nThe quasi-Hamiltonian is defined as in Eq. (20) of Joynt et al. This function returns -i H_q, in order to deal with real matrices instead of imaginary ones.\n\nArguments\n\nHamiltonian is an array of matrices, each of them representing the Hamiltonian\n\nwith a different configuration of the noise.\n\nγ is the switching rate of the noise.\n\n\n\n\n\n"
},

{
    "location": "#QuasiHamiltonianRTN.Vnoise",
    "page": "QuasiHamiltonianRTN.jl",
    "title": "QuasiHamiltonianRTN.Vnoise",
    "category": "function",
    "text": "Vnoise(n, γ)\n\nReturns the generator V for the transition probabilities between the various noise configurations for n fluctuators, with switching rate gamma.\n\nIt is defined as\n\nmathbfV = sum_i=1^N mathbfV_gamma^(i)\n\nwith mathbfV_gamma^(i) = mathbfI_2^otimes i-1 otimes mathbfV_gamma otimes mathbfI_2^otimes N-i.\n\n\n\n\n\n"
},

{
    "location": "#QuasiHamiltonianRTN.commutator",
    "page": "QuasiHamiltonianRTN.jl",
    "title": "QuasiHamiltonianRTN.commutator",
    "category": "function",
    "text": "commutator(A, B)\n\nReturns the commutator [A,B] = A * B - B * A between two matrices\n\n\n\n\n\n"
},

{
    "location": "#QuasiHamiltonianRTN.density_operator",
    "page": "QuasiHamiltonianRTN.jl",
    "title": "QuasiHamiltonianRTN.density_operator",
    "category": "function",
    "text": "density_operator(r)\n\nReturns the density operator ρ corresponding to the Bloch vector r.\n\nvec must be a vector of N^2-1 real values. We follow the formula ρ = frac1N (I + sqrtN mathbfr cdot mathbfλ)\n\n\n\n\n\n"
},

{
    "location": "#QuasiHamiltonianRTN.random_density_matrix",
    "page": "QuasiHamiltonianRTN.jl",
    "title": "QuasiHamiltonianRTN.random_density_matrix",
    "category": "function",
    "text": "random_density(n)\n\nReturns a random N  N density matrix, distributed according to the Hilbert-Schmidt measure.\n\nNOTE: No guarantee about the quality of the algorithm For more details see V Al Osipov et al 2010 J. Phys. A: Math. Theor. 43 055302\n\n\n\n\n\n"
},

{
    "location": "#QuasiHamiltonianRTN.structure_constants",
    "page": "QuasiHamiltonianRTN.jl",
    "title": "QuasiHamiltonianRTN.structure_constants",
    "category": "function",
    "text": "structure_constants(λ)\n\nReturns the totally antisymmetric tensor i f_ijk of the structure constants of the list of generators λ, defined as f_ijk = - 2 i textTr(λ_iλ_jλ_k).\n\nNotice that the true definition has a factor -i, which we avoid to store only a real value.\n\n\n\n\n\n"
},

{
    "location": "#QuasiHamiltonianRTN.unitary_op",
    "page": "QuasiHamiltonianRTN.jl",
    "title": "QuasiHamiltonianRTN.unitary_op",
    "category": "function",
    "text": "unitary_op(H, t)\n\nReturns the unitary evolution operator U = exp-i H t.\n\n\n\n\n\n"
},

{
    "location": "#Reference-1",
    "page": "QuasiHamiltonianRTN.jl",
    "title": "Reference",
    "category": "section",
    "text": "QuasiHamiltonianRTN.evolution\nQuasiHamiltonianRTN.igen\nQuasiHamiltonianRTN.bloch_vector\nQuasiHamiltonianRTN.quasiHamiltonian\nQuasiHamiltonianRTN.Vnoise\nQuasiHamiltonianRTN.commutator\nQuasiHamiltonianRTN.density_operator\nQuasiHamiltonianRTN.random_density_matrix\nQuasiHamiltonianRTN.structure_constants\nQuasiHamiltonianRTN.unitary_op"
},

{
    "location": "method/#",
    "page": "Details on the method",
    "title": "Details on the method",
    "category": "page",
    "text": ""
},

{
    "location": "method/#Details-on-the-method-1",
    "page": "Details on the method",
    "title": "Details on the method",
    "category": "section",
    "text": "Joynt et al. have proposed an exact method of solving the dynamics of a quantum system coupled to a classical environment modeled as a Markovian stochastic process, and particularly effective for RTN. The method allows for analytical results only for a single qubit, while it requires numerical matrix diagonalization for higher dimensions. The strengths of this method are that it gives exact results up to machine-precision, and it avoids fluctuations typical of Montecarlo simulations: the drawback, however, is the exponential complexity in terms of the number of noise fluctuators.Suppose to have a N_q-dimensional quantum state, described at time t by the N_q times N_q density matrix rho(t), and a classical system made of N_c states, representing the possible values of the noise. For example, if the classical noise is a single fluctuator, N_c = 2; if it consists of N independent RTN sources, N_c = 2^N.We start at t=0 with rho(0) and the classical probability distribution mathbfP(0), describing the initial state of the stochastic process associated to the classical noise. The Hamiltonian of the quantum system is Hg(t), i.e., a function of the stochastic process describing the noise. At every time instant, to every particular configuration of the noise, which we label with the index c in 1 ldots N_c, corresponds a particular form of the Hamiltonian H_c.Since we assume a Markovian classical environment, the transitions between different states c c is described by a master equationbeginequation\nfrac d mathbfP(t) dt = mathbfV mathbfP(t)\nendequationwhere the element mathbfV_cc dictates the transition rates between the states c and c of the environment. Notice that mathbfV is time-independent if the stochastic process describing the environment is homogeneous, as is the case in this work. In the case of a single RTN with a switching rate mu, we have that N_c = 2 andbeginequation\nlabeleqVrtnMU\nmathbfV_mu =beginpmatrix\n-mu  mu \nmu  -mu\nendpmatrix\nendequationFor a collection of N independent fluctuators, the matrix mathbfV becomesbeginequation\nlabeleqVrtnN\n  mathbfV = sum_i=1^N mathbfV_mu^(i) quad mathbfV_mu^(i) = mathbfI_2^otimes i-1 otimes mathbfV_mu otimes mathbfI_2^otimes N-i\nendequationWe \\eqref{eq:VrtnMU} \\eqref{eq:VrtnN} need to represent the density matrix as a vector, and we do so by employing the generalized Bloch vector mathbfn(t), a vector of dimension N_q^2-1, with real componentsbeginequation\n  newcommandTrtextTr\n  n_i(t) = fracsqrtN_q2 Tr lambda_i rho(t)\nendequationwhere lambda_j are the generators of SU(N_q), and they are N_qtimes N_q matrices chosen to satisfybeginequation\nlabeleqSUNgenerators\nTr lambda_j =0quad lambda_j^dagger=lambda_jquad Tr (lambda_jlambda_k) = 2 delta_jk\nendequationWe can go back to the density matrix rho(t) from the Bloch vector mathbfn(t) with the equationbeginequation\nlabeleqSUNmixedState\nrho(t)=frac1N_qleftmathbfI_N_q+sqrtN_qsum_j=1^N_q^2-1 n_j(t)lambda_jright\nendequationwhere mathbfI_N_q denotes the identity in the Hilbert space of the quantum system.The action of a unitary operator U onto the density matrix rho(t) is translated into the multiplication of the Bloch vector mathbfn(t) by a transfer matrix T defined asbeginequation\n  mathbfT_ij = frac 12 Trlambda_i U lambda_j U^dagger\nendequationConsider now a short time interval Delta t in which the environment is in a fixed state c; during Delta t, the unitary evolution is generated by the Hamiltonian H_c: U_c(Delta t)=exp -iH_c Delta t. The corresponding transfer matrix mathbfT_c is generated by the matrixbeginequationlabeleqquasihamiltonian_generators\n  G_c = - lim_Delta t to 0 fracmathbfT_c - mathbfI_N_qDelta t = - frac 1 2 sum_ij=1^N_q Tr left(lambda_i lambda_j H_c right)\nendequationJoynt et al. introduced the quasi-Hamiltonian matrixbeginequation\nlabeleqquasiHamDivided\nnewcommandHqmathbfH_q\nHq= - mathbfVotimesmathbfI_N_q^2-1 + bigoplus_i=c^N_c G_c\nendequationwhere the second term is a direct sum of all the generators defined in \\eqref{eq:quasihamiltonian_generators}, and showed that the dynamics of the system, averaged all the possible realizations of the stochastic process describing the noise (as defined in Eq. \\eqref{eq:map}), is given bybeginequationlabeleqhqfinalequation\n  newcommandket1left 1 rightrangle\n  mathbfn(t) = leftlanglemathbf1bigexp(- Hq t)bigp_0rightranglecdotmathbfn(0)\nendequationIn Eq. \\eqref{eq:hqfinalequation}, ketp_0 and ketmathbf1 are vectors belonging to the space of the classical configurations: ketmathbf1 is a vector with all components set to 1, while ketp_0 equiv mathbfP(0) is the initial probability distribution of the configurations of the noise. In the case at hand, where we assume stationary noise, all the configurations are equally probable and sobeginequation\n  ketp_0 = frac 1 N_c ketmathbf1\nendequationThe expression langlemathbf1Ap_0rangle where A is a N_c (N_q^2-1) times N_c (N_q^2-1) matrix, denotes a partial inner product in the space of classical configurations: the result is a (N_q^2 -1)  times (N_q^2 -1) matrix acting on the Bloch vector of the quantum system.With the definitions chosen above, which are slightly different from the original paper, the matrix Hq is real, which allows us to avoid the use of complex numbers at all."
},

]}
