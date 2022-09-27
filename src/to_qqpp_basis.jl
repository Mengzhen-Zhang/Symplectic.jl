export toQQPPBasis

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
