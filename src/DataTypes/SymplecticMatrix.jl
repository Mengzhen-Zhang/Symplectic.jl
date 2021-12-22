# Define the Symplectic matrix data type
struct SymplecticMatrix{T} <: AbstractMatrix{T}
    matrix:: AbstractMatrix{T}
    modes:: Int
    function SymplecticMatrix{T}(matrix::AbstractMatrix) where T
        if !(T<:Real) return error("not real") end
        if !(LinearAlgebraUtilities.isSymplectic(matrix)) return error("not symplectic") end
        return new{T}(AbstractMatrix{T}(matrix), size(matrix)[1] ÷ 2)
    end
end

SymplecticMatrix(matrix::AbstractMatrix) = SymplecticMatrix{typeof(matrix).parameters[1]}(matrix)
SymplecticMatrix(S::SymplecticMatrix) = S
Omega(n::Int)::SymplecticMatrix = LinearAlgebraUtilities.Omega(n)
Omega(S::SymplecticMatrix)::SymplecticMatrix = LinearAlgebraUtilities.Omega(S.matrix)
Ω = Omega

isGeneric(S::SymplecticMatrix)::Bool = LinearAlgebraUtilities.isGeneric(S.matrix)


# Override functions in Base
Base.Matrix(S::SymplecticMatrix)::AbstractMatrix = S.matrix
Base.convert(::Type{SymplecticMatrix}, matrix::AbstractMatrix) = SymplecticMatrix(matrix)
Base.convert(::Type{SymplecticMatrix}, S::SymplecticMatrix) = S
Base.convert(::Type{AbstractMatrix}, S::SymplecticMatrix) = Matrix(S)
Base.size(S::SymplecticMatrix, dim::Integer) = size(S.matrix, dim)
Base.size(S::SymplecticMatrix) = size(S.matrix)
## Import !!!
function Base.getindex(S::SymplecticMatrix{T}, i::Integer, j::Integer)::T where {T} 
    return S.matrix[i, j]
end
Base.transpose(S::SymplecticMatrix)::SymplecticMatrix = Base.transpose( S.matrix )
Base.adjoint(S::SymplecticMatrix)::SymplecticMatrix = Base.adjoint( S.matrix )
# Unary Minus and Unary Plus
Base.:-(S::SymplecticMatrix)::SymplecticMatrix = Base.:-( S.matrix )
Base.:+(S::SymplecticMatrix)::SymplecticMatrix = Base.:+( S.matrix )
Base.:*(S1::SymplecticMatrix, S2::SymplecticMatrix, Ss::Vararg{SymplecticMatrix})::SymplecticMatrix = *(S1.matrix, S2.matrix,  [s.matrix for s in Ss]... )
Base.inv(S::SymplecticMatrix)::SymplecticMatrix = - Ω(S) * S' * Ω(S) 
Base.:^(S::SymplecticMatrix, p::Int)::SymplecticMatrix = p < 0 ? inv(S).matrix ^ (-p)  : S.matrix^p
⊗ = Base.kron

# Direct Sum
dsum(M1::AbstractMatrix, M2::AbstractMatrix)::AbstractMatrix = LinearAlgebraUtilities.dsum(M1, M2)
dsum(S1::SymplecticMatrix, S2::SymplecticMatrix, Ss::Vararg{SymplecticMatrix})::SymplecticMatrix = dsum(S1.matrix, S2.matrix,  [s.matrix for s in Ss]... )
⊕ = dsum
