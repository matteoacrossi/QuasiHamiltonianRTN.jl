# Details on the method

Joynt et al. have proposed an exact method of solving the dynamics of a quantum system coupled to a classical environment modeled as a Markovian stochastic process, and particularly effective for RTN. The method allows for analytical results only for a single qubit, while it requires numerical matrix diagonalization for higher dimensions. The strengths of this method are that it gives exact results up to machine-precision, and it avoids fluctuations typical of Montecarlo simulations: the drawback, however, is the exponential complexity in terms of the number of noise fluctuators.

Suppose to have a $N_q$-dimensional quantum state, described at time $t$ by the $N_q \times N_q$ density matrix $\rho(t)$, and a classical system made of $N_c$ states, representing the possible values of the noise. For example, if the classical noise is a single fluctuator, $N_c = 2$; if it consists of $N$ independent RTN sources, $N_c = 2^N$.

We start at $t=0$ with $\rho(0)$ and the classical probability distribution $\mathbf{P}(0)$, describing the initial state of the stochastic process associated to the classical noise. The Hamiltonian of the quantum system is $H[g(t)]$, i.e., a function of the stochastic process describing the noise. At every time instant, to every particular configuration of the noise, which we label with the index $c \in \{1, \ldots, N_c\}$, corresponds a particular form of the Hamiltonian $H_c$.

Since we assume a Markovian classical environment, the transitions between different states $c, c'$ is described by a master equation

```math
\begin{equation}
\frac {d \mathbf{P}(t)} {dt} = \mathbf{V} \mathbf{P}(t),
\end{equation}
```

where the element $\mathbf{V}_{c,c'}$ dictates the transition rates between the states $c$ and $c'$ of the environment. Notice that $\mathbf{V}$ is time-independent if the stochastic process describing the environment is homogeneous, as is the case in this work. In the case of a single RTN with a switching rate $\mu$, we have that $N_c = 2$ and
```math
\begin{equation}
\label{eq:VrtnMU}
\mathbf{V}_\mu =\begin{pmatrix}
-\mu & \mu \\
\mu & -\mu
\end{pmatrix}.
\end{equation}
```

For a collection of $N$ independent fluctuators, the matrix $\mathbf{V}$ becomes
```math
\begin{equation}
\label{eq:VrtnN}
  \mathbf{V} = \sum_{i=1}^N \mathbf{V}_\mu^{(i)}, \quad \mathbf{V}_\mu^{(i)} = \mathbf{I}_2^{\otimes i-1} \otimes \mathbf{V}_\mu \otimes \mathbf{I}_2^{\otimes N-i}.
\end{equation}
```
We \eqref{eq:VrtnMU} \eqref{eq:VrtnN} need to represent the density matrix as a vector, and we do so by employing the generalized Bloch vector $\mathbf{n}(t)$, a vector of dimension $N_q^2-1$, with real components
```math
\begin{equation}
  \newcommand{\Tr}{\text{Tr}}
  n_i(t) = \frac{\sqrt{N_q}}{2} \Tr \lambda_i \rho(t),
\end{equation}
```
where $\lambda_j$ are the generators of $SU(N_q)$, and they are $N_q\times N_q$ matrices chosen to satisfy
```math
\begin{equation}
\label{eq:SUNgenerators}
\Tr \lambda_j =0,\quad \lambda_j^\dagger=\lambda_j,\quad \Tr (\lambda_j\lambda_k) = 2 \delta_{jk}.
\end{equation}
```
We can go back to the density matrix $\rho(t)$ from the Bloch vector $\mathbf{n}(t)$ with the equation
```math
\begin{equation}
\label{eq:SUNmixedState}
\rho(t)=\frac{1}{N_q}\left[\mathbf{I}_{N_q}+\sqrt{N_q}\sum_{j=1}^{N_q^2-1} n_j(t)\lambda_j\right],
\end{equation}
```
where $\mathbf{I}_{N_q}$ denotes the identity in the Hilbert space of the quantum system.

The action of a unitary operator $U$ onto the density matrix $\rho(t)$ is translated into the multiplication of the Bloch vector $\mathbf{n}(t)$
by a transfer matrix $T$ defined as
```math
\begin{equation}
  \mathbf{T}_{ij} = \frac 12 \Tr[\lambda_i U \lambda_j U^\dagger].
\end{equation}
```
Consider now a short time interval $\Delta t$ in which the environment is in a fixed state $c$; during $\Delta t$, the unitary evolution is generated by the Hamiltonian $H_c$: $U_c(\Delta t)=\exp [-iH_c \Delta t]$. The corresponding transfer matrix $\mathbf{T}_c$ is generated by the matrix
```math
\begin{equation}\label{eq:quasihamiltonian_generators}
  G_c = - \lim_{\Delta t \to 0} \frac{\mathbf{T}_c - \mathbf{I}_{N_q}}{\Delta t} = - \frac 1 2 \sum_{i,j=1}^{N_q} \Tr \left([\lambda_i, \lambda_j] H_c \right).
\end{equation}
```
[Joynt *et al.*](https://dx.doi.org/10.1142/S0217979211100990) introduced the *quasi-Hamiltonian* matrix
```math
\begin{equation}
\label{eq:quasiHamDivided}
\newcommand{\Hq}{\mathbf{H}_q}
\Hq= - \mathbf{V}\otimes\mathbf{I}_{N_q^2-1} + \bigoplus_{i=c}^{N_c} G_c,
\end{equation}
```
where the second term is a direct sum of all the generators defined in \eqref{eq:quasihamiltonian_generators},
and showed that the dynamics of the system, averaged all the possible realizations of the stochastic process describing the noise (as defined in Eq. \eqref{eq:map}), is given by
```math
\begin{equation}\label{eq:hqfinalequation}
  \newcommand{\ket}[1]{\left| #1 \right\rangle}
  \mathbf{n}(t) = \left\langle{\mathbf{1}}\big|{\exp(- \Hq t)}\big|{p_0}\right\rangle\cdot\mathbf{n}(0).
\end{equation}
```

In Eq. \eqref{eq:hqfinalequation}, $\ket{p_0}$ and $\ket{\mathbf{1}}$ are vectors belonging to the space of the classical configurations: $\ket{\mathbf{1}}$ is a vector with all components set to $1$, while $\ket{p_0} \equiv \mathbf{P}(0)$ is the initial probability distribution of the configurations of the noise. In the case at hand, where we assume stationary noise, all the configurations are equally probable and so
```math
\begin{equation}
  \ket{p_0} = \frac 1 {N_c} \ket{\mathbf{1}}.
\end{equation}
```
The expression $\langle\mathbf{1}|{A}|{p_0}\rangle$ where $A$ is a $N_c (N_q^2-1) \times N_c (N_q^2-1)$ matrix,
denotes a partial inner product in the space of classical configurations: the result
is a $(N_q^2 -1)  \times (N_q^2 -1)$ matrix acting on the Bloch vector of the quantum system.

With the definitions chosen above, which are slightly different from the original
paper, the matrix $\Hq$ is real, which allows us to avoid the use of complex numbers at all.
