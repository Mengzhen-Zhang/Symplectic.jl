export dsum, ⊕

"""
    dsum(M1, M2, M3...)

Return the matrix direct sum of the arguments.
"""
dsum(M1::AbstractMatrix, M2::AbstractMatrix) = cat(M1, M2; dims=(1, 2))
dsum(M1::AbstractMatrix, M2::AbstractMatrix, Ms::Vararg{AbstractMatrix}) = dsum(dsum(M1, M2), mapreduce(identity, dsum, Ms))
const ⊕ = dsum

