using LinearAlgebra

# Define the AbstractSymplecticMatrix data type
abstract type AbstractSymplecticMatrix{T} <: AbstractArray{T, 2} end
Base.size(A::AbstractSymplecticMatrix) = error("size not implemented")
Base.getindex(A::AbstractSymplecticMatrix, i::Int) = error("getindex not implemented")
Base.getindex(A::AbstractSymplecticMatrix, I::Vararg{Int, N}) = error("getindex not implemented")
Base.setindex!(A::AbstractSymplecticMatrix, v, i::Int) = error("setindex! not implemented")
Base.setindex!(A::AbstractSymplecticMatrix, v, I::Vararg{Int, N}) = error("setindex! not implemented")
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
Base.getindex(S::SymplecticMatrix, I::Vararg{Int, N}) = getindex(S.matrix, I)
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
Base.convert(::Type{SymplecticMatrix}, M::AbstractArray{T, 2}) = SymplecticMatrix(M)
Base.convert(::Type{AbstractArray}, S::SymplecticMatrix) = Matrix(S)


###########################################################
### TODO: Define SymplectifForm as a separate data type ###
###########################################################

struct SymplecticForm{T<:Number} 
    λ::T
end

const Ω = SymplecticForm(1)

(Ω::UniformScaling)(n::Integer) = cat(fill([zero(Ω.λ) Ω.λ; -Ω.λ zero(Ω.λ)], n)...; dims=(1,2))
Base.eltype(::Type{SymplecticForm{T}}) where {T} = T
Base.ndims(J::SymplecticForm) = 2
Base.has_offset_axes(::UniformScaling) = false
function Base.getindex(J::SymplecticForm, i::Integer, j::Integer)
    if i==j-1
        return J.λ
    elseif i==j+1
        return -J.λ
    else
        return zero(J.λ)
    end
end

Base.getindex(J::SymplecticForm, n::Integer, m::AbstractRange{<:Integer}) = -getindex(x, m, n)
function Base.getindex(J::SymplecticForm{T}, n::AbstractRange{<:Integer}, m::Integer) where T
    v = zeros(T, length(n))
    @inbounds for (i, ii) in enumerate(n)
        if ii == m - 1
            v[i] = -J.λ
        elseif ii == m + 1
            v[i] = J.λ
        end
    end
    return v
end

function Base.getindex(J::UniformScaling{T}, n::AbstractRange{<:Integer}, m::AbstractRange{<:Integer}) where T
    A = zeros(T, length(n), length(m))
    @inbounds for (j, jj) in enumerate(m), (i, ii) in enumerate(n)
        if ii == jj - 1
            v[i] = -J.λ
        elseif ii == jj + 1
            v[i] = J.λ
        end
    end
end

function Base.show(io::IO, ::MIME"text/plain", J::SymplecticForm)
    s = "$(J.λ)"
    if occursin(r"\w+\s*[\+\-]\s*\w+", s)
        s="($s)"
    end
    print(io, typeof(J), "\n$s*Ω")
end
Base.copy(J::SymplecticForm) = SymplecticForm(J.λ)

Base.convert(::Type{SymplecticForm{T}}, J::SymplecticForm) where{T} = SymplecticForm(convert(T, J.λ))

Base.transpose(J::SymplecticForm) = J
Base.adjoint(J::SymplecticForm) = transpose(J)

Base,one(::Type{SymplecticForm{T}}) where {T} = SymplecticForm(one(T))
Base.one(J::SymplecticForm{T}) where {T} = one(SymplecticForm{T})
Base.oneunit(::Type{SymplecticForm{T}}) where {T} = SymplecticForm(oneunit(T))
Base.oneunit(J::SymplecticForm{T}) where {T} = oneunit(SymplecticForm{T})
Base.zero(::Type{SymplecticForm{T}}) where {T} = SymplecticForm(zero(T))
Base.zero(J::SymplecticForm{T}) where {T} = zero(SymplecticForm{T})

LinearAlgebra.isdiag(::SymplecticForm) = false
LinearAlgebra.istriu(::SymplecticForm) = false
LinearAlgebra.istril(::SymplecticForm) = false
LinearAlgebra.issymmetric(::SymplecticForm) = false
LinearAlgebra.ishermitian(::SymplecticForm) = false
LinearAlgebra.isposdef(::SymplecticForm) = false

(Base.:+)(J::SymplecticForm) = SymplecticForm(+J.λ)
(Base.:+)(J1::SymplecticForm, J2::SymplecticForm) = SymplecticForm(J1.λ + J2.λ)
(Base.:+)(J::SymplecticForm, A::AbstractMatrix) = A + J
function (Base.:+)(A::AbstractMatrix, J::SymplecticForm)
    checksquare(A)
    B = convert(AbstractMatrix{Base._return_type(+, Tuple{eltype(A), eltype(J)})}, A)
    for i in 1:size(B)[1]-1
        @inbounds B[i,i+1] += J.λ
        @inbounds B[i+1, i] -= J.λ
    end
    return B
end


(Base.:-)(J::SymplecticForm) = SymplecticForm(-J.λ)
(Base.:-)(J1::SymplecticForm, J2::SymplecticForm) = SymplecticForm(J1.λ - J2.λ)
(Base.:-)(A::AbstractMatrix, J::SymplecticForm) = A + (-J)
function (Base.:-)(J::SymplecticForm, A::AbstractMatrix)
    checksquare(A)
    B = convert(AbstractMatrix{Base._return_type(+, Tuple{eltype(A), eltype(J)})}, -A)
    for i in 1:size(B)[1]-1
        @inbounds B[i,i+1] += J.λ
        @inbounds B[i+1, i] -= J.λ
    end
    return B
end

Base.inv(J::SymplecticForm) = SymplecticForm(-inv(J.λ))

LinearAlgebra.opnorm(J::SymplecticForm, p::Real=2) = opnorm(J.λ, p)

LinearAlgebra.pinv(J::SymplecticForm) = iszero(J.λ) ? SymplecticForm(zero(inv(J.λ))) : SymplecticForm(-inv(J.λ))

function det(J::SymplecticForm{T}) where T
    if isone(J.λ)
        one(T)
    elseif iszero(J.λ)
        zero(T)
    else
        throw(ArgumentError("Determinant of SymplecticForm is only well-defined when λ=0 or 1."))
    end
end

LinearAlgebra.tr(J::SymplecticForm{T}) where T = zero(T)

Base.:*(J1::SymplecticForm, J2::SymplecticForm) = SymplecticForm(J1.λ * J2.λ)
function (Base.:*)(A::AbstractMatrix, J::SymplecticForm)
    if size(A)[2] % 2 != 0
        throw(DimensionMismatch("SymplecticForm is always of even dimension."))
    end
    B = zeros(eltype(J), size(A)[1], size(A)[2])
    for j in 1:size(B)[2]
        @inbounds B[:, j] = j % 2 == 0 ? J.λ * A[:, j-1] : (-J.λ) * A[:, j+1]
    end
    return B
end
function (Base.:*)(J::SymplecticForm, A::AbstractMatrix)
    if size(A)[1] % 2 != 0
        throw(DimensionMismatch("SymplecticForm is always of even dimension."))
    end
    B = zeros(eltype(J), size(A)[1], size(A)[2])
    for i in 1:size(B)[1]
        @inbounds B[i, :] = i % 2 == 0 ? J.λ * A[i-1, :] : (-J.λ) * A[i+1, :]
    end
    return B
end
(Base.:*)(v::AbstractVector, J::SymplecticForm) = reshape(v, 1, length(v)) * J
(Base.:*)(J::SymplecticForm, v::AbstractVector) = J * reshape(v, length(v), 1)
(Base.:*)(x::Number, J::SymplecticForm) = SymplecticForm(x * J.λ)
(Base.:*)(J::SymplecticForm, x::Number) = SymplecticForm(J.λ * x)

(Base.:/)(J1::SymplecticForm, J2::SymplecticForm) = J2.λ == 0 ? throw(SingularException(1)) : UniformScaling(J1.λ / J2.λ)
(Base.:/)(J::SymplecticForm, A::AbstractMatrix) = J * inv(A)
(Base.:/)(A::AbstractMatrix, J::SymplecticForm) = J.λ == 0 ? throw(SingularException(1)) : A * inv(J)
(Base.:/)(v::Abstractvector, J::SymplecticForm) = reshape(v, length(v), 1) * inv(J)
(Base.:/)(J::SymplecticForm, x::Number) = SymplecticForm(J.λ / x)

(Base.:\)(J1::SymplecticForm, J2::SymplecticForm) = J1.λ == 0 ? throw(SingularException(1)) : SymplecticForm(J1.λ \ J2.λ)
(Base.:\)(J::SymplecticForm, A::AbstractVecOrMat) = J.λ == 0 ? throw(SingularException(1)) : A / J
(Base.:\)(A::AbstractMatrix, J::SymplecticForm) = J / A
(Base.:\)(F::Factorization, J::SymplecticForm) = F \ J(size(F, 1))
(Base.:\)(x::Number, J::SymplecticForm) = SymplecticForm(x \ J.λ)








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

