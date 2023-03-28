using LinearAlgebra
using ChainRulesCore

export inverseSymplecticCayleyTransform, invcayley
export symplecticCayleyTransform, cayley
export nonSymplecticity
export feedforward
export teleport
export interferenceBasedSequence, intfseq
export circulator
export toQPQPBasis
export toQQPPBasis
export preIwasawaFactorization, preiwa
export spectrum
export dilate
export adaptiveMeasurement

export dsum, ⊕
export phaseShifting, ps
export amplifier, amp
export extend
export beamSplitter, bs
export randomPhaseShifiting, randps
export localSymplecticMatrix, localsymp
export randomLocalSymplecticMatrix, randlocalsymp
export ⊗
export squeezedVacuum, channel

"""
    inverseSymplecticCayleyTransform(S::AbstractMatrix)

return `-Ω*((I-S)^-1 - I/2)`.
"""
inverseSymplecticCayleyTransform(S::AbstractMatrix) = -Ω * ((I - S)^-1 - I / 2)
const invcayley = inverseSymplecticCayleyTransform

"""
    symplecticCayleyTransform(M::AbstractMatrix)

Return `I - (Ω * M + I / 2)^-1`.
"""
symplecticCayleyTransform(M::AbstractMatrix) = I - (Ω * M + I / 2)^-1
const cayley = symplecticCayleyTransform

"""
    nonSymplecticity(A)

Return the matrix norm of ``AΩA^T-Ω``.

`A` is symplectic when return zero. Work for `SymplecticForm` and `UniformScaling` as well.
"""
nonSymplecticity(J::Union{SymplecticForm,UniformScaling}) = 0
function nonSymplecticity(A::AbstractMatrix)
       n = size(A, 1)
       M = transpose(A) * Ω * A - Ω * I(n)
       tr(transpose(M) * M)
end

"""
    feedforward(S::AbstractMatrix, inModes::Vector, outModes::Vector)

Caculate Eq. (3) in *PHYSICAL REVIEW LETTERS 120, 020502 (2018)*, with modes not contained in
`inModes` squeezed along their Q-quadratures and modes not contained in `outModes` measured
along their P-quadratures.
"""
function feedforward(S::AbstractMatrix, inModes::Vector, outModes::Vector)
       if length(inModes) == 0
              return S
       end
       allModes = [1:(size(S, 1)÷2);]
       In = vcat([[2 * i - 1, 2 * i] for i in inModes]...)
       Out = vcat([[2 * i - 1, 2 * i] for i in outModes]...)
       ancModes = [i for i in allModes if !(i in inModes)]
       idlModes = [i for i in allModes if !(i in outModes)]
       Usq = 2 * ancModes .- 1
       Hm = 2 * idlModes .- 1
       SOutUsq = S[Out, Usq]
       SHmUsq = S[Hm, Usq]
       return - SOutUsq * Matrix(SHmUsq)^-1
end


"""
    teleportation(S::AbstractMatrix, inModes::Vector, outModes::Vector)

Caculate Eq. (4) in *PHYSICAL REVIEW LETTERS 120, 020502 (2018)*, with modes not contained in
`inModes` squeezed along their Q-quadratures and modes not contained in `outModes` measured
along their P-quadratures.
"""
function teleport(S::AbstractMatrix, inModes::Vector, outModes::Vector)
       if length(inModes) == 0
              return S
       end
       allModes = [1:(size(S, 1)÷2);]
       In = vcat([[2 * i - 1, 2 * i] for i in inModes]...)
       Out = vcat([[2 * i - 1, 2 * i] for i in outModes]...)
       ancModes = [i for i in allModes if !(i in inModes)]
       idlModes = [i for i in allModes if !(i in outModes)]
       Usq = 2 * ancModes .-1
       Hm = 2 * idlModes .- 1
       SOutIn = S[Out, In]
       SHmIn = S[Hm, In]
       SOutUsq = S[Out, Usq]
       SHmUsq = S[Hm, Usq]
       return SOutIn - SOutUsq * Matrix(SHmUsq)^-1 * SHmIn
end

# ↑(m, i) = 2 * m - (i % 2)
# ↓(m, i) = 2 * m - 1 + (i % 2)

"""
    Lkmij(u, v, m, i, j)

Compute Eq.(14) in *npj Quantum Information 8 (1), 1-10*.
"""
function Lkmij(u, v, m, i, j)
       if norm(v[2*m-1:2*m]) ≈ 0
           return i == j ? 1 : 0
       else
           up_i = 2 * m - (i % 2)
           up_j = 2 * m - (j % 2)
           down_i = 2 * m - 1 + (i % 2)
           down_j = 2 * m - 1 + (j % 2)
           return (-1)^(j + 1) * v[up_i] * u[down_j] / (v[2*m]^2 + v[2*m-1]^2) + (-1)^i * v[down_i] * u[up_j] / (u[2*m]^2 + u[2*m-1]^2)
       end
end

"""
    localSymplecticTransformation(u::Vector, v::Vector, lambda::Real=1)

Return a local symplectic transformation transforming quadrature u to quadrature λΩv (the conjugate of quadrature v)
"""
function localSymplecticTransformation(u::Vector, v::Vector; lambda::Real=1)
       (length(u) == length(v) && length(u) % 2 == 0 && length(v) % 2 == 0) || throw(ArgumentError("invalid input"))
       n = length(u) ÷ 2
       if n == 1
              return [Lkmij(u, lambda * v, 1, i, j) for i in 1:2, j in 1:2]
       else
              return ⊕([[Lkmij(u, lambda * v, m, i, j) for i in 1:2, j in 1:2] for m in 1:n]...)
       end
end
const lstf = localSymplecticTransformation


"""
    decoupleSequence(S4, S3, S2 S1; m::Int=1, lambda::Real=1)

Return a sequence S₄ L₃ S₃ L₂ S₂ L₁ S₁ decoupling the m-th mode, where the
local operations are determined up to the scaling factor λ.
"""
function decoupleSequence(S4, S3, S2, S1; m::Int=1, lambda::Real=1)
       L1 = lstf(S1[:, 2*m-1], S2[2*m-1, :]; lambda=lambda)
       L3 = lstf(S3[:, 2*m-1], S4[2*m-1, :]; lambda=lambda^-1)
       T1 = S2 * L1 * S1
       T2 = S4 * L3 * S3
       L2 = ⊕([(i == m ? Ω(2) : lstf(T1[2*i-1:2*i, 2*m], T2[2*m, 2*i-1:2*i])) for i in 1:(size(T1, 1)÷2)]...)
       return [S4, L3, S3, L2, S2, L1, S1]
end

"""
   getSequence(Ss, modeToDecouple, lambda::Real=1)

Return a vetor of symplectic matrices with the modes specificed by `modeToDecouple` decoupled.

`Ss` is a vector of given symplectic matrices. If `modeToDecouple=m`, then resulting sequence
decouples the first `m` modes.
"""
function getSequence(Ss, modeToDecouple::Int; lambda::Real=1)
       if length(Ss) == 4
              return decoupleSequence(Ss...; m=modeToDecouple, lambda=lambda)
       else
              step = length(Ss) ÷ 4
              Rs = [getSequence(Ss[i:i+step-1], modeToDecouple - 1) for i in 1:step:length(Ss)]

              Ls = decoupleSequence([*(RLst...) for RLst in Rs]...; m=modeToDecouple, lambda=lambda)
              return vcat([(i % 2 == 0) ? [Ls[i]] : Rs[(i+1)÷2] for i in 1:7]...)
       end
end


"""
    interferenceBasedSequence(S; T=I(4), lambda::Real=1)

Return the sequence as in the right-hand side of Eq. (1) of *npj Quantum Information 8 (1), 1-10*.

Get local Operations L16, L15, ..., L1 for the interference-basded sequence
LR S L15 S ... S L1 S = ST, with a given target Symplectic Matrix ST, if S is single matrix. S can
also be a vector of symplectic matrices.
"""
function interferenceBasedSequence(S::AbstractMatrix; T=I(4), lambda::Real=1)
       return interferenceBasedSequence([S for i in 1:4^(size(T,1)÷2)]; T=T, lambda=lambda)
end
function interferenceBasedSequence(Ss; T=I(4), lambda::Real=1)
       n = (size(Ss[1], 1)÷2) - (size(T, 1)÷2)
       nT = size(T, 1)÷2
       Ss[end] = Ss[end] * (inv(T) ⊕ I(2 * n))
       Ss = getSequence(Ss, nT; lambda=lambda)
       R = (*(Ss...)[1:2*nT, 1:2*nT])^-1 ⊕ I(2 * n)
       Ss[end] = Ss[end] * (T ⊕ I(2 * n))
       return vcat([R], Ss)
end
const intfseq = interferenceBasedSequence

"""
    circulator(perm...)

Return a symplectic matrix permutating modes.

Input can be a vector.
"""
circulator(perm...) = circulator(perm)
circulator(perm::Vector) =
       hvcat(length(perm), [i == perm[j] ? I(2) : zeros(2, 2) for i in 1:length(perm), j in 1:length(perm)]...)



"""
    toQPQPBasis(mat::AbstractMatrix)

Return a matrix represented in the Q1, P1, Q2, P2, ... basis.

`mat` is assumed to be a  matrix in the Q1, Q2, ...; P1, P2, ... basis.
"""
function toQPQPBasis(mat::AbstractMatrix)
    m, n = size(mat)
    ord_m = vcat([[i, i + (m ÷ 2)] for i in 1:(m ÷ 2)]...)
    if length(ord_m) < m
        push!(ord_m, m)
    end
    ord_n = vcat([[i, i + (n ÷ 2)] for i in 1:(n ÷ 2)]...)
    if length(ord_n) < n
        push!(ord_n, n)
    end

    return mat[ord_m, ord_n]
end


"""
    toQQPPBasis(mat::AbstractMatrix)

Return a symplectic matrix represented in the Q1, Q2,...; P1, P2, ... basis.

`mat` is assumed to be a  matrix in the Q1, P1, Q2, P2, ... basis.
"""
function toQQPPBasis(mat::AbstractMatrix)
    m, n = size(mat)
    ord_m = [[1:2:m;]..., [2:2:m;]...]
    ord_n = [[1:2:n;]..., [2:2:m;]...]
    return mat[ord_m, ord_n]
end


"""
    preIwasawaFactorization(S::AbstractMatrix)

Return vector of three matrices representing the pre-Iwasawa-factorization of a symplectic matrix.

The product of the output is equal to `S`.
"""
function preIwasawaFactorization(S::AbstractMatrix)
       n = size(S, 1)÷2
       S = toQQPPBasis(S)
       A, B = S[1:n, 1:n], S[1:n, 1+n:2*n]
       C, D = S[1+n:2*n, 1:n], S[1+n:2*n, 1+n:2*n]
       AAtBBt = A * transpose(A) + B * transpose(B)
       P = (C * transpose(A) + D * transpose(B)) * AAtBBt^-1
       L = AAtBBt^(1 / 2)
       X = AAtBBt^(-1 / 2) * A
       Y = AAtBBt^(-1 / 2) * B
       return map(toQPQPBasis, [
              [I(n) zeros(n, n); P I(n)],
              [L zeros(n, n); zeros(n, n) transpose(L)^-1],
              [X Y; -Y X]
       ])
end
const preiwa = preIwasawaFactorization

"""
    spectrum(M::AbstractMatrix)

Return the symplectic spectrum of `M`, in the form of [vector of values, matrix `S`].

Sᵀ M S = toQPQPBasis( Diagonal([λ..., λ...]) )
"""
function spectrum(M::AbstractMatrix)
       (eltype(M) <: Real && isposdef(M)) || throw(ArgumentError("only implemented for real-valued positive-definite matrix"))
       n = size(M, 1) ÷ 2
       K = toQQPPBasis(Ω * M)
       eig = eigen(K; sortby=(λ -> (abs(λ), imag(λ))))
       spec = [abs(imag(v)) for v in eig.values]
       es = [(v = eig.vectors[:, 2*i-1]; real(v)) for i in 1:n]
       fs = [(v = eig.vectors[:, 2*i-1]; -imag(v)) for i in 1:n]
       for i in 1:n
              λ = es[i] ⋅ (toQQPPBasis(Ω(2 * n)) * fs[i])
              es[i] /= sqrt(λ)
              fs[i] /= sqrt(λ)
       end
       S = hcat(es..., fs...)
       return [spec, toQPQPBasis(S)]
end


"""
    dsum(M1, M2, M3...)

Return the matrix direct sum of the arguments.
"""
dsum(M1::AbstractMatrix, M2::AbstractMatrix) = cat(M1, M2; dims=(1, 2))
dsum(M1::AbstractMatrix, M2::AbstractMatrix, Ms::Vararg{AbstractMatrix}) = dsum(dsum(M1, M2), mapreduce(identity, dsum, Ms))
const ⊕ = dsum


"""
    phaseShifting(θ...)

Return a (multi-)mode phase-shifting symplectic matrix.

Size of the returned matrix is specified by the length of the input.
"""
phaseShifting(θ::Real) = [cos(θ) sin(θ); -sin(θ) cos(θ)]
phaseShifting(θs...) = dsum(map(phaseShifting, θs)...)
const ps = phaseShifting


"""
    amplifier(G::Real[, mode1, mode2[, modes]])

Return a two-mode-squeezing symplectic matrix, acting on `mode1` and `mode2` of
a system containg `modes` of modes.

Return default to return a 4×4 matrix. `modes` default to `max(mode1, mode2)`.
"""
amplifier(G::Real) = [sqrt(G) 0 sqrt(G-1) 0; 0 sqrt(G) 0 -sqrt(G-1); sqrt(G-1) 0 sqrt(G) 0; 0 -sqrt(G-1) 0 sqrt(G)]
# Generate a BeamSplitter between two modes in a multi-mode System
function amplifier(G::Real, mode1, mode2, modes)
       mode1, mode2 = sort([mode1, mode2])
       lst = [[[2 * (i + 2) - 1, 2 * (i + 2)] for i in 1:mode1-1]...; [1, 2]; [[2 * (i + mode1 + 1) - 1, 2 * (i + mode1 + 1)] for i in 1:mode2-mode1-1]...; [3, 4]; [[2 * (i + mode2) - 1, 2 * (i + mode2)] for i in 1:modes-mode2]...]
       return dsum(amplifier(G), I(2 * modes - 4))[lst, lst]
end
amplifier(G::Real, mode1, mode2) = amplifier(G, mode1, mode2, max(mode1, mode2))
const amp = amplifier

"""
    extend(m::AbstractMatrix, pos::AbstractVector{Int}, modes::Int)

Return a direct sum, `A`, of `m` and an identity matrix, with the `pos[i]`th mode of `A` equal to
the `i`th mode of the input matrix `m`.
"""
function extend(m::AbstractMatrix, pos::AbstractVector{Int}, modes::Int)
    lst = zeros(Int64, modes)
    for (i, p) in enumerate(pos)
        lst[p] = i
    end
    j = size(m, 1) + 1
    for i in 1:modes
        if lst[i] == 0
            lst[i] = j
            j += 1
        end
    end
    return dsum(m, I(modes - size(m, 1)))[lst, lst]
end


"""
    beamSplitter(angle::Real[, mode1, mode2[, modes]])

Return a beam-splitter symplectic matrix, acting on `mode1` and `mode2` of
a system containg `modes` of modes.

Return default to return a 4×4 matrix. `modes` default to `max(mode1, mode2)`.
"""
beamSplitter(angle::Real) = kron(phaseShifting(angle), I(2))
# Generate a BeamSplitter between two modes in a multi-mode System
function beamSplitter(angle::Real, mode1, mode2, modes)
       mode1, mode2 = sort([mode1, mode2])
       lst = [[[2 * (i + 2) - 1, 2 * (i + 2)] for i in 1:mode1-1]...; [1, 2]; [[2 * (i + mode1 + 1) - 1, 2 * (i + mode1 + 1)] for i in 1:mode2-mode1-1]...; [3, 4]; [[2 * (i + mode2) - 1, 2 * (i + mode2)] for i in 1:modes-mode2]...]
       return dsum(beamSplitter(angle), I(2 * modes - 4))[lst, lst]
end
beamSplitter(angle::Real, mode1, mode2) = beamSplitter(angle, mode1, mode2, max(mode1, mode2))
const bs = beamSplitter


"""
    randomPhaseShifiting(n::Int)

Return an n-mode random phase-shifting symplectic matrix.
"""
randomPhaseShifiting(n::Int) =
       n == 1 ? phaseShifting(2π * rand()) : phaseShifting((2π * rand(n))...)
const randps = randomPhaseShifiting

localSymplecticMatrix(θ::Real, l::Real, ϕ::Real) =
    phaseShifting(θ) * [l 0; 0 1/l] * phaseShifting(ϕ)

"""
    localSymplecticMatrix(θs::Union{Real,Vector}, ls::Union{Real,Vector}, ϕs::Union{Real,Vector})

Return a direct sum of 2×2 `localSymplecticMatrix`s according to Euler decompostion.

Return a 2×2 matrix when the input arguments are real numbers.
"""
localSymplecticMatrix(θs::Vector, ls::Vector, ϕs::Vector) =
    dsum(map(args -> localSymplecticMatrix(args...), zip(θs, ls, ϕs))...)
const localsymp = localSymplecticMatrix


"""
    randomLocalSymplecticMatrix(n::Int, bound::Real=0.2)

Return a random direct sum of 2×2 `localSymplecticMatrix`s.

`bound` bounds the degree of squeezing.
"""
randomLocalSymplecticMatrix(n::Int; bound::Real=0.2) =
       n == 1 ? localSymplecticMatrix(2π * rand(), bound * (rand() - 1 / 2) + 1.0, 2π * rand()) :
       localSymplecticMatrix(2π * rand(n), bound * (rand(n) .- 1 / 2) .+ 1.0, 2π * rand(n))
const randlocalsymp = randomLocalSymplecticMatrix

⊗ = Base.kron

function dilate(S::AbstractMatrix)
    M = invcayley(S)
    Mₐ = (M - transpose(M)) / 2
    chol = ignore_derivatives() do
           skewchol(2 * Ω * Mₐ * Ω) 
        end
    R = ignore_derivatives() do 
        transpose(chol.R[:, invperm(chol.p)])
    end 
    L = - Ω * transpose(R) * Ω
    S1 = I - L * inv(Ω * M + I / 2) * R
    return [S (I - S)*R; L*(I - S) S1]
end

function adaptiveMeasurement(F::AbstractMatrix, outModes::Vector, modes::Integer)
    rowF, colF = size(F)
    E = transpose(F) * Ω
    G = - transpose(F) * Ω * F / 2
    order = vcat([[i, i + colF] for i in 1:colF]...)
    F1 = [F   zeros(rowF, colF)][:, order]
    E1 = [zeros(colF, rowF); -E][order, :]
    G1 = toQPQPBasis([I(colF)  zeros(colF, colF); G  I(colF)])
    A = [I(rowF)     F1;
         E1          G1]
    other = [m for m in 1:modes if m ∉ outModes]
    order = vcat(outModes, other)
    perm = ChainRulesCore.ignore_derivatives() do 
        vcat([[2*p-1, 2*p] for p in invperm(order)]...) 
    end
    return A[perm, perm]
end

function channel(S::AbstractMatrix, Venv::AbstractMatrix, inModes::Vector, outModes::Vector)
    n = size(S, 1) ÷ 2
    inSys = vcat([[2*m-1, 2*m] for m in inModes]...)
    outSys = vcat([[2*m-1, 2*m] for m in outModes]...)
    inEnv = vcat([[2*m-1, 2*m] for m in 1:n if m ∉ inModes]...)
    T = S[outSys, inSys]
    N = S[outSys, inEnv] * Venv * transpose(S[outSys, inEnv])
    return (T, N)
end

function squeezedVacuum(ξs)
    return diagm(vcat([[exp(-ξ * log(10)), exp(ξ * log(10))] for ξ in ξs]...))
end
