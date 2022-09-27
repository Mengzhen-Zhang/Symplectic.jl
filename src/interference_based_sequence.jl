using LinearAlgebra: I, norm

export interferenceBasedSequence, ibs

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
const ibs = interferenceBasedSequence
