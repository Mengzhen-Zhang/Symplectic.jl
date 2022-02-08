# Utilities
using LinearAlgebra: checksquare, opnorm

include("SymplecticForm.jl")

# calculate the 'nonSymplecticity".
nonSymplecticity(A::AbstractMatrix) = opnorm(A'*Ω*A - Ω) / opnorm(A)

# check if a matrix is symplectic, returning its non-symplecticity
function checkSymplectic(A::AbstractMatrix)
    A <: Real || throw(ArgumentError("matrix is not real-valued"))
    checksquare(A) % 2 == 0 || throw(DimensionMismatch("matrix is not even-dimensional"))
    return nonSymplecticity(A)
end
function checkSymplecticMatrix(A...)
    nonSymplecticities = []
    for a in A
        a <: Real || throw(ArgumentError("matrix is not real-valued"))
        checksquare(a) % 2 == 0 || throw(DimensionMismatch("matrix is not even-dimensional"))
        push!(nonSymplecticities, nonSymplecticity(a))
    end
    return nonSymplecticities
end







    


using LinearAlgebra

# Define the AbstractSymplecticMatrix data type
abstract type AbstractSymplecticMatrix{T} <: AbstractArray{T, 2} end
Base.size(A::AbstractSymplecticMatrix) = error("size not implemented")
Base.getindex(A::AbstractSymplecticMatrix, i::Int) = error("getindex not implemented")
Base.getindex(A::AbstractSymplecticMatrix, I::Vararg{Int, N}) where {N} = error("getindex not implemented")
Base.setindex!(A::AbstractSymplecticMatrix, v, i::Int) = error("setindex! not implemented")
Base.setindex!(A::AbstractSymplecticMatrix, v, I::Vararg{Int, N}) where {N} = error("setindex! not implemented")
# get a Q-quadrature from columns
getQCol(A::AbstractSymplecticMatrix, i::Int) = error("getQCol not implemented")
# get a P-quadrature from columns
getPCol(A::AbstractSymplecticMatrix, i::Int) = error("getPCol not implemented")
# get a Q-quadrature from rows
getQRow(A::AbstractSymplecticMatrix, i::Int) = error("getQCol not implemented")
# get a P-quadrature from rows
getPRow(A::AbstractSymplecticMatrix, i::Int) = error("getPCol not implemented")

Base.transpose(S::AbstractSymplecticMatrix) = error("transpose not implemented")
Base.adjoint(S::AbstractSymplecticMatrix) = error("adjoin not implemented")

Base.:-(S::AbstractSymplecticMatrix) = error("- not implemented")
Base.:+(S::AbstractSymplecticMatrix) = error("+ not implemented")
Base.:*(S1::AbstractSymplecticMatrix, S2::AbstractSymplecticMatrix, Ss::Vararg{AbstractSymplecticMatrix}) = error("* not implemented")
Base.inv(S::AbstractSymplecticMatrix) = error("inv not implemented")
Base.:^(S::AbstractSymplecticMatrix, p::Int) = error("^ not implemented")


# Define the SymplecticMatrix data type
struct SymplecticMatrix{T} <: AbstractSymplecticMatrix{T}
    matrix:: AbstractArray{T, 2}
    function SymplecticMatrix{T}(matrix:: AbstractArray{T, 2}) where T
        if !(T<:Real) return error("not real") end
        if !(LinearAlgebraUtilities.isSymplectic(matrix)) return error("not symplectic") end
        return new{T}(AbstractMatrix{T}(matrix))
    end
end

Base.size(S::SymplecticMatrix) = size(S.matrix)
Base.getindex(S::SymplecticMatrix, i::Int) = getindex(S.matrix, i)
Base.getindex(S::SymplecticMatrix, I::Vararg{Int, N}) where {N} = getindex(S.matrix, I)
getQCol(S::SymplecticMatrix, i::Int) = S.matrix[:, 2*i-1]
getPCol(S::SymplecticMatrix, i::Int) = S.matrix[:, 2*i]
getQRow(S::SymplecticMatrix, i::Int) = S.matrix[2*i-1, :]
getPRow(S::SymplecticMatrix, i::Int) = S.matrix[2*i, :]

SymplecticMatrix(M::AbstractArray{T, 2}) where T = SymplecticMatrix{T}(M)
SymplecticMatrix(S::SymplecticMatrix) = S

Base.transpose(S::SymplecticMatrix)::SymplecticMatrix = transpose( S.matrix )
Base.adjoint(S::SymplecticMatrix)::SymplecticMatrix = adjoint( S.matrix )

Base.:-(S::SymplecticMatrix)::SymplecticMatrix = Base.:-( S.matrix )
Base.:+(S::SymplecticMatrix)::SymplecticMatrix = Base.:+( S.matrix )
Base.:*(S1::SymplecticMatrix, S2::SymplecticMatrix, Ss::Vararg{SymplecticMatrix})::SymplecticMatrix = *(S1.matrix, S2.matrix,  [s.matrix for s in Ss]... )
Base.inv(S::SymplecticMatrix)::SymplecticMatrix = - Ω(S) * S' * Ω(S) 
Base.:^(S::SymplecticMatrix, p::Int)::SymplecticMatrix = p < 0 ? inv(S).matrix ^ (-p)  : S.matrix^p

# Unary Minus and Unary Plus
⊗ = Base.kron

# Direct Sum
dsum(M1::AbstractMatrix, M2::AbstractMatrix)::AbstractMatrix = LinearAlgebraUtilities.dsum(M1, M2)
dsum(S1::SymplecticMatrix, S2::SymplecticMatrix, Ss::Vararg{SymplecticMatrix})::SymplecticMatrix = dsum(S1.matrix, S2.matrix,  [s.matrix for s in Ss]... )
⊕ = dsum

Base.Matrix(S::SymplecticMatrix)::Array{T, 2} = S.matrix
Base.convert(::Type{SymplecticMatrix}, S::SymplecticMatrix) = S
Base.convert(::Type{SymplecticMatrix}, M::AbstractArray{T, 2}) where T = SymplecticMatrix(M)
Base.convert(::Type{AbstractArray}, S::SymplecticMatrix) = Matrix(S)


# Omega(n::Int)::SymplecticMatrix = LinearAlgebraUtilities.Omega(n)
# Omega(S::SymplecticMatrix)::SymplecticMatrix = LinearAlgebraUtilities.Omega(S.matrix)
# Ω = Omega

isGeneric(S::SymplecticMatrix)::Bool = LinearAlgebraUtilities.isGeneric(S.matrix)


# Override functions in Base
# Base.Matrix(S::SymplecticMatrix)::AbstractMatrix = S.matrix
# Base.convert(::Type{SymplecticMatrix}, matrix::AbstractMatrix) = SymplecticMatrix(matrix)
# Base.convert(::Type{SymplecticMatrix}, S::SymplecticMatrix) = S
# Base.convert(::Type{AbstractMatrix}, S::SymplecticMatrix) = Matrix(S)
# Base.size(S::SymplecticMatrix, dim::Integer) = size(S.matrix, dim)
# Base.size(S::SymplecticMatrix) = size(S.matrix)
## Import !!!
# function Base.getindex(S::SymplecticMatrix{T}, i::Integer, j::Integer)::T where {T} 
#     return S.matrix[i, j]
# end

