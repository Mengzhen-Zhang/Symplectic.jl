# depend on SymplecticForm.jl

using LinearAlgebra: tr, UniformScaling

export nonSymplecticity

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
