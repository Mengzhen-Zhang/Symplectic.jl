module Symplectic

import LinearAlgebra

include("DataTypes/SymplecticForm.jl")

export SymplecticForm, Ω, nonSymplecticity, checkSymplectic, modes,
       qQuadrature, pQuadrature, ⊗, dsum, ⊕

# Utilities
using LinearAlgebra: checksquare, opnorm

# calculate the 'nonSymplecticity".
nonSymplecticity(A::AbstractMatrix) = opnorm(A'*Ω*A - Ω) / opnorm(Ω)

# check if a matrix is symplectic, returning its non-symplecticity
function checkSymplectic(A::AbstractMatrix)
       A <: Real || throw(ArgumentError("matrix is not real-valued"))
       checksquare(A) % 2 == 0 || throw(DimensionMismatch("matrix is not even-dimensional"))
       return nonSymplecticity(A)
end
function checkSymplectic(A...)
       nonSymplecticities = []
       for a in A
              a <: Real || throw(ArgumentError("matrix is not real-valued"))
              checksquare(a) % 2 == 0 || throw(DimensionMismatch("matrix is not even-dimensional"))
              push!(nonSymplecticities, nonSymplecticity(a))
       end
       return nonSymplecticities
end

modes(A::AbstractMatrix) = checkSymplectic(A) && size(A, 1) ÷ 2

# Get quadtrature (column vector) from a symplectic matrix
qQuadrature(A::AbstractMatrix, n::Int) = 
       n ≤ modes(A) ? A[:, 2*n-1] : throw(ArgumentError("index is out of bound"))
function qQuadrature(A::AbstractMatrix{T}, n::AbstractRange{<:Integer}) where T
       length(n) ≤ modes(A) || throw(ArgumentError("index is out of bound"))
       @inbounds v = [A[:, 2*i-1] for i in n]
       return hcat(v...)
end

pQuadrature(A::AbstractMatrix, n::Int) = 
       n ≤ modes(A) ? A[:, 2*n] : throw(ArgumentError("index is out of bound"))
function pQuadrature(A::AbstractMatrix{T}, n::AbstractRange{<:Integer}) where T
       length(n) ≤ modes(A) || throw(ArgumentError("index is out of bound"))
       @inbounds v = [A[:, 2*i-1] for i in n]
       return hcat(v...)
end

# const ⊗ = Base.kron

# dsum(M1::AbstractMatrix, M2::AbstractMatrix)::AbstractMatrix = cat(M1, M2; dims=(1,2))
# dsum(M1::AbstractMatrix, M2::AbstractMatrix, Ms::Vararg{AbstractMatrix})::AbstractMatrix = dsum(dsum(M1, M2),  mapreduce(identity, dsum, Ms ))
# ⊕ = dsum







    
#     include("Utilities/LinearAlgebraUtilities.jl")
#     import .LinearAlgebraUtilities

#     include("DataTypes/SymplecticMatrix.jl")

#     include("DataTypes/MatrixSequence.jl")



#     # Data type: SymplecticMatrix
#     export SymplecticMatrix,
#            Omega, Ω,
#            isGeneric,
#            ⊗,
#            dsum, ⊕

#     # Data type: MatrixSequence
#     export MatrixSequence,
#            sliceInner,
#            slice,
#            setInterspersion,
#            setSequence,
#            replaceSequence

#     # Utility module for linear algebra
#     export LinearAlgebraUtilities

end