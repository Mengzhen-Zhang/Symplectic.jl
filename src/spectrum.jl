# depend on to_qqpp_basis.jl
# depend on to_qpqp_basis.jl

using LinearAlgebra: isposdef, eigen

export spectrum

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
