export toQPQPBasis

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
