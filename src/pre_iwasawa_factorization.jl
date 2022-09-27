# depend on to_qqpp_basis.jl
# depend on to_qpqp_basis.jl

export preIwasawaFactorization, pif

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
              [L zeros(n, n); zeros(n, n) L^-1],
              [X Y; -Y X]
       ])
end
const pif = preIwasawaFactorization
