# Extended Methods
import Base: eltype, ndims, has_offset_axes, getindex, show, copy,
             convert, transpose, adjoint, one, oneunit, zero, +,
             -, inv, *, /, \, ^, ==, isequal, isapprox, hcat, vcat,
             hvcat, Matrix, Array, muladd, conj, real, imag

import LinearAlgebra: isdiag, istriu, istril, issymmetric, ishermitian,
                    isposdef, opnorm, pinv, det, tr, mul!, rmul!, lmul!,
                    rdiv!, ldiv!, copyto!, cond, promote_to_arrays_,
                    promote_to_array_type, Diagonal, dot, _diag_or_value

# Utilities
using LinearAlgebra: checksquare, Factorization, promote_to_arrays, I, 
                    UniformScaling
using Base: require_one_based_indexing

using ChainRulesCore

export SymplecticForm, Ω

struct SymplecticForm{T<:Number} 
    λ::T
end

const Ω = SymplecticForm(1)

function (Ω::SymplecticForm)(n::Integer) 
    B = zeros(Int, n, n)
    @inbounds for i in 1:n
        if i % 2 == 0
            B[i, i-1] = -1
        elseif i + 1 ≤ n
            B[i, i+1] = 1
        end
    end
    return B
end

eltype(::Type{SymplecticForm{T}}) where {T} = T

ndims(J::SymplecticForm) = 2

has_offset_axes(::SymplecticForm) = false

function getindex(J::SymplecticForm, i::Integer, j::Integer)
    if i % 2 == 0 && j == i - 1 
        return -J.λ 
    elseif i % 2 == 1 && j == i + 1 
        return J.λ 
    else 
        return zero(J.λ)
    end
end
getindex(J::SymplecticForm, n::Integer, m::AbstractRange{<:Integer}) = -getindex(x, m, n)
function getindex(J::SymplecticForm{T}, n::AbstractRange{<:Integer}, m::Integer) where T
    v = zeros(T, length(n))
    @inbounds for (i, ii) in enumerate(n)
        if ii % 2 == 0 && m == ii - 1
            v[i] = -J.λ
        elseif ii % 2 == 1 && m == ii + 1
            v[i] = J.λ
        else
            v[i] = zero(J.λ)
        end
    end
    return v
end
function getindex(J::SymplecticForm{T}, n::AbstractRange{<:Integer}, m::AbstractRange{<:Integer}) where T
    A = zeros(T, length(n), length(m))
    @inbounds for (j, jj) in enumerate(m), (i, ii) in enumerate(n)
        if ii == jj - 1 && ii % 2 == 0
            A[i, j] = -J.λ
        elseif ii == jj + 1 && ii % 2 == 1
            A[i ,j] = J.λ
        else
            A[i, j] = zero(J.λ)
        end
    end
    return A
end

function show(io::IO, ::MIME"text/plain", J::SymplecticForm)
    s = "$(J.λ)"
    if occursin(r"\w+\s*[\+\-]\s*\w+", s)
        s="($s)"
    end
    print(io, typeof(J), "\n$s*Ω")
end

copy(J::SymplecticForm) = SymplecticForm(J.λ)

convert(::Type{SymplecticForm{T}}, J::SymplecticForm) where {T} = SymplecticForm(convert(T, J.λ))

conj(J::SymplecticForm) = SymplecticForm(conj(J.λ))

real(J::SymplecticForm) = SymplecticForm(real(J.λ))

imag(J::SymplecticForm) = SymplecticForm(imag(J.λ))

transpose(J::SymplecticForm) = SymplecticForm(- J.λ)

adjoint(J::SymplecticForm) = SymplecticForm(- conj(J.λ))

one(::Type{SymplecticForm{T}}) where {T} = UniformScaling(one(T))
one(J::SymplecticForm{T}) where {T} = one(SymplecticForm{T})

oneunit(::Type{SymplecticForm{T}}) where {T} = UniformScaling(oneunit(T))
oneunit(J::SymplecticForm{T}) where {T} = oneunit(SymplecticForm{T})

zero(::Type{SymplecticForm{T}}) where {T} = SymplecticForm(zero(T))
zero(J::SymplecticForm{T}) where {T} = zero(SymplecticForm{T})

isdiag(::SymplecticForm) = false

istriu(::SymplecticForm) = false

istril(::SymplecticForm) = false

issymmetric(::SymplecticForm) = false

ishermitian(::SymplecticForm) = false

isposdef(::SymplecticForm) = false

(+)(J::SymplecticForm, x::Number)           = zero(J.λ) + x
(+)(x::Number, J::SymplecticForm)           = x + zero(J.λ)
(+)(J::SymplecticForm)                      = SymplecticForm(+J.λ)
(+)(J::SymplecticForm, B::BitArray{2})      = J + Array(B)
(+)(B::BitArray{2}, J::SymplecticForm)      = Array(B) + J
(+)(J1::SymplecticForm, J2::SymplecticForm) = SymplecticForm(J1.λ + J2.λ)
(+)(J::SymplecticForm, A::AbstractMatrix)   = A + J
function (+)(A::AbstractMatrix, J::SymplecticForm)
    checksquare(A)
    B = convert(AbstractMatrix{Base._return_type(+, Tuple{eltype(A), eltype(J)})}, A)
    @inbounds for i in 1:size(B, 1)
        if i % 2 == 0
            B[i, i - 1] -= J.λ
        elseif i + 1 ≤ size(B, 2)
            B[i, i + 1] += J.λ
        end
    end
    return B
end

(-)(J::SymplecticForm, x::Number)           = zero(J.λ) - x
(-)(x::Number, J::SymplecticForm)           = x - zero(J.λ)
(-)(J::SymplecticForm)                      = SymplecticForm(-J.λ)
(-)(J1::SymplecticForm, J2::SymplecticForm) = SymplecticForm(J1.λ - J2.λ)
(-)(B::BitArray{2}, J::SymplecticForm)      = Array(B) - J
(-)(J::SymplecticForm, B::BitArray{2})      = J - Array(B)
(-)(A::AbstractMatrix, J::SymplecticForm)   = A + (-J)
function (-)(J::SymplecticForm, A::AbstractMatrix)
    checksquare(A)
    B = convert(AbstractMatrix{Base._return_type(+, Tuple{eltype(A), eltype(J)})}, -A)
    @inbounds for i in 1:size(B, 1)
        if i % 2 == 0
            B[i, i - 1] -= J.λ
        elseif i + 1 ≤ size(B, 2)
            B[i, i + 1] += J.λ
        end
    end
    return B
end

inv(J::SymplecticForm) = SymplecticForm(-inv(J.λ))

opnorm(J::SymplecticForm, p::Real=2) = opnorm(J.λ, p)

pinv(J::SymplecticForm) = iszero(J.λ) ? SymplecticForm(zero(inv(J.λ))) : SymplecticForm(-inv(J.λ))

function det(J::SymplecticForm{T}) where T
    if isone(J.λ)
        one(T)
    elseif iszero(J.λ)
        zero(T)
    else
        throw(ArgumentError("Determinant of SymplecticForm is only well-defined when λ=0 or 1."))
    end
end

tr(J::SymplecticForm{T}) where T = zero(T)

(*)(J::SymplecticForm, U::UniformScaling) = SymplecticForm(J.λ * U.λ)
(*)(U::UniformScaling, J::SymplecticForm) = SymplecticForm(U.λ * J.λ)
(*)(J1::SymplecticForm, J2::SymplecticForm) = UniformScaling(- J1.λ * J2.λ)
(*)(J::SymplecticForm, B::BitArray{2}) = J * Array(B)
(*)(B::BitArray{2}, J::SymplecticForm) = Array(B) * J
function (*)(A::AbstractMatrix, J::SymplecticForm)
    B = zeros(Base._return_type(*, Tuple{eltype(A), eltype(J)}), size(A))
    ignore_derivatives() do 
        @inbounds for j in 1:size(B, 2)
            if j % 2 == 0
                B[:, j] = A[:, j-1] * J.λ
            elseif j + 1 ≤ size(B, 2)
                B[:, j] = - A[:, j+1] * J.λ
            end
        end
    end
    return B
end
function (*)(J::SymplecticForm, A::AbstractMatrix)
    B = zeros(Base._return_type(*, Tuple{eltype(J), eltype(A)}), size(A))
    @inbounds for i in 1:size(B, 1)
        if i % 2 == 0
            B[i, :] = - J.λ * A[i - 1, :]
        elseif i + 1 ≤ size(B, 1)
            B[i, :] = J.λ * A[i + 1, :]
        end
    end
    return B
end
function (*)(v::AbstractVector, J::SymplecticForm)
    w = zeros(Base._return_type(*, Tuple{eltype(v), eltype(J)}), size(v))
    @inbounds for i in 1:size(v, 1)
        if i % 2 == 0
            w[i] = v[i-1] * J.λ
        elseif i + 1 ≤ size(v, 1)
            w[i] = -v[i+1] * J.λ
        end
    end
    return w
end
function (*)(J::SymplecticForm, v::AbstractVector)
    w = zeros(Base._return_type(*, Tuple{eltype(J), eltype(v)}), size(v))
    @inbounds for i in 1:size(v, 1)
        if i % 2 == 0
            w[i] = -J.λ * v[i-1]
        elseif i + 1 ≤ size(v, 1)
            w[i] = J.λ * v[i+1]
        end
    end
    return w
end
(*)(x::Number, J::SymplecticForm) = SymplecticForm(x * J.λ)
(*)(J::SymplecticForm, x::Number) = SymplecticForm(J.λ * x)

(/)(J1::SymplecticForm, J2::SymplecticForm) = J2.λ == 0 ? throw(SingularException(1)) : UniformScaling(J1.λ / J2.λ)
(/)(J::SymplecticForm, A::AbstractMatrix) = J * inv(A)
(/)(A::AbstractMatrix, J::SymplecticForm) = J.λ == 0 || size(A, 2) % 2 == 1 ? throw(SingularException(1)) : A * inv(J)
(/)(v::AbstractVector, J::SymplecticForm) = J.λ == 0 || size(v, 1) % 2 == 1 ? throw(SingularException(1)) : v * inv(J)
(/)(J::SymplecticForm, x::Number) = SymplecticForm(J.λ / x)

(\)(J1::SymplecticForm, J2::SymplecticForm) = J1.λ == 0 ? throw(SingularException(1)) : UniformScaling(J1.λ \ J2.λ)
(\)(J::SymplecticForm, A::AbstractVecOrMat) = J.λ == 0 || size(A, 1) % 2 == 1 ? throw(SingularException(1)) : inv(J) * A
(\)(A::AbstractMatrix, J::SymplecticForm) = inv(A) * J
(\)(F::Factorization, J::SymplecticForm) = F \ J(size(F, 1))
(\)(x::Number, J::SymplecticForm) = SymplecticForm(x \ J.λ)

@inline mul!(C::AbstractMatrix, A::AbstractMatrix, J::SymplecticForm, alpha::Number, beta::Number) =
    mul!(C, I(size(A, 1)),  A*J, alpha, beta)
@inline mul!(C::AbstractMatrix, J::SymplecticForm, B::AbstractVecOrMat, alpha::Number, beta::Number) = 
    mul!(C, J*B, I(size(B, 2)), alpha, beta)
function mul!(out::AbstractMatrix{T}, a::Number, J::SymplecticForm, alpha::Number, beta::Number) where {T}
    checksquare(out)
    if iszero(beta)
        fill!(out, zero(T))
    elseif !isone(beta)
        rmul!(out, beta)
    end
    s = convert(T, a*J.λ*alpha)
    if !iszero(s)
        @inbounds for i in 1:size(out, 1)
            if i % 2 == 0
                out[i, i-1] += s
            elseif i + 1 ≤ size(out, 2)
                out[i, i+1] -= s
            end
        end
    end
    return out
end
@inline mul!(out::AbstractMatrix, J::SymplecticForm, b::Number, alpha::Number, beta::Number) = 
    mul!(out, J.λ, SymplecticForm(b), alpha, beta)
rmul!(A::AbstractMatrix, J::SymplecticForm) = mul!(A, A, J)
lmul!(J::SymplecticForm, B::AbstractVecOrMat) = mul!(B, J, B)
rdiv!(A::AbstractMatrix, J::SymplecticForm) = size(A, 2) % 2 == 1 ? throw(SingularException(1)) : mul!(A, A, inv(J))
ldiv!(J::SymplecticForm, B::AbstractVecOrMat) = size(B, 1) % 2 == 1 ? throw(SingularException(1)) : mul!(B, inv(J), B)
ldiv!(Y::AbstractVecOrMat, J::SymplecticForm, B::AbstractVecOrMat) = (Y .= J \ B)

Broadcast.broadcasted(::typeof(*), x::Number, J::SymplecticForm) = SymplecticForm(x * J.λ)
Broadcast.broadcasted(::typeof(*), J::SymplecticForm, x::Number) = SymplecticForm(J.λ * x)

Broadcast.broadcasted(::typeof(/), J::SymplecticForm, x::Number) = SymplecticForm(J.λ / x)

Broadcast.broadcasted(::typeof(\), x::Number, J::SymplecticForm) = SymplecticForm(x \ J.λ)

function (^)(J::SymplecticForm, x::Int)
    if x >= 0
        return x % 2 == 0 ? UniformScaling((-J.λ)^(x ÷ 2)) : (-J.λ)^(x ÷ 2) * J
    else
        return inv(J) ^ (- x)
    end
end
Broadcast.broadcasted(::typeof(^), J::SymplecticForm, x::Number) = SymplecticForm(J.λ ^ x)

(==)(J1::SymplecticForm, J2::SymplecticForm) = (J1.λ == J2.λ)
(==)(J::SymplecticForm, A::AbstractMatrix) = A == J
function ==(A::AbstractMatrix, J::SymplecticForm)
    try
        checksquare(A)
    catch e
        return false
    end
    if iszero(J.λ) return iszero(A) end
    return A == J(size(A, 1))
end

isequal(A::AbstractMatrix, J::SymplecticForm) = false
isequal(J::SymplecticForm, A::AbstractMatrix) = false

function isapprox(J1::SymplecticForm{T}, J2::SymplecticForm{S}; 
                atol::Real=0, rtol::Real=Base.rtoldefault(T, S, atol), nans::Bool=false) where {T<:Number, S<:Number}
    return isapprox(J1.λ, J2.λ, rtol=rtol, atol=atol, nans=nans)
end
function isapprox(J::SymplecticForm, A::AbstractMatrix;
                atol::Real=0,
                rtol::Real=Base.rtoldefault(promote_leaf_eltypes(A), eltype(J), atol),
                nans::Bool=false, norm::Function=norm
                )
    n = checksquare(A)
    normJ = norm === opnorm             ? abs(J.λ) : 
            norm === LinearAlgebra.norm ? abs(J.λ) * sqrt(n) :
                                          norm(cat(fill([zero(J.λ) J.λ; -J.λ zero(J.λ)], n)...; dims=(1,2)))
    return norm(A - J) <= max(atol, rtol*max(norm(A), normJ))
end
isapprox(A::AbstractMatrix, J::SymplecticForm; kwargs...) = isapprox(J, A; kwargs...)

function copyto!(A::AbstractMatrix, J::SymplecticForm)
    Base.require_one_based_indexing(A)
    fill!(A, 0)
    λ = J.λ
    s = min(size(A, 1), size(A, 2))
    @inbounds for i = 1:s
        if i % 2 == 0
            A[i, i - 1] = - λ
        elseif i + 1 ≤ size(A, 2)
            A[i ,i + 1] = λ
        end
    end
    return A
end

function cond(J::SymplecticForm{T}) where T
    onereal = inv(one(real(J.λ)))
    return J.λ ≠ zero(T) ? onereal : oftype(onereal, Inf)
end

promote_to_arrays_(n::Int, ::Type{Matrix}, J::SymplecticForm{T}) where{T} = copyto!(Matrix{T}(undef, n, n), J)
promote_to_array_type(J::Tuple{Vararg{Union{AbstractVecOrMat, UniformScaling, SymplecticForm}}}) = Matrix

for (f, dim, name) in ((:hcat, 1, "rows"), (:vcat, 2, "cols"))
    @eval begin
        function $f(A::Union{AbstractVecOrMat, UniformScaling, SymplecticForm}...)
            n = -1
            for a in A
                if !isa(a, Union{SymplecticForm, UniformScaling})
                    na = size(a, $dim)
                    n >= 0 && n != na &&
                        throw(DimensionMismatch(string("number of ", $name,
                            " of each array must match (got ", n, " and ", na, ")")))
                    n = na
                end
            end
            n == -1 && throw(ArgumentError($("$f of only SymplecticForm or UniformScaling objects cannot determine the matrix size")))
            return $f(promote_to_arrays(fill(n, length(A)), 1, promote_to_array_type(A), A...)...)
        end
    end
end

function hvcat(rows::Tuple{Vararg{Int}}, A::Union{AbstractVecOrMat,UniformScaling, SymplecticForm}...)
    require_one_based_indexing(A...)
    nr = length(rows)
    sum(rows) == length(A) || throw(ArgumentError("mismatch between row sizes and number of arguments"))
    n = fill(-1, length(A))
    needcols = false # whether we also need to infer some sizes from the column count
    j = 0
    for i = 1:nr # infer UniformScaling sizes from row counts, if possible:
        ni = -1 # number of rows in this block-row, -1 indicates unknown
        for k = 1:rows[i]
            if !isa(A[j+k], Union{UniformScaling, SymplecticForm})
                na = size(A[j+k], 1)
                ni >= 0 && ni != na &&
                    throw(DimensionMismatch("mismatch in number of rows"))
                ni = na
            end
        end
        if ni >= 0
            for k = 1:rows[i]
                n[j+k] = ni
            end
        else # row consisted only of UniformScaling and/or SymplecticForm objects
            needcols = true
        end
        j += rows[i]
    end
    if needcols # some sizes still unknown, try to infer from column count
        nc = -1
        j = 0
        for i = 1:nr
            nci = 0
            rows[i] > 0 && n[j+1] == -1 && (j += rows[i]; continue)
            for k = 1:rows[i]
                nci += isa(A[j+k], Union{UniformScaling, SymplecticForm}) ? n[j+k] : size(A[j+k], 2)
            end
            nc >= 0 && nc != nci && throw(DimensionMismatch("mismatch in number of columns"))
            nc = nci
            j += rows[i]
        end
        nc == -1 && throw(ArgumentError("sizes of UniformScalings could not be inferred"))
        j = 0
        for i = 1:nr
            if rows[i] > 0 && n[j+1] == -1 # this row consists entirely of UniformScalings
                nci = nc ÷ rows[i]
                nci * rows[i] != nc && throw(DimensionMismatch("indivisible UniformScaling sizes"))
                for k = 1:rows[i]
                    n[j+k] = nci
                end
            end
            j += rows[i]
        end
    end
    return hvcat(rows, promote_to_arrays(n,1, promote_to_array_type(A), A...)...)
end

Matrix(J::SymplecticForm, dims::Dims{2}) = Matrix{eltype(J)}(J, dims)
function Matrix{T}(J::SymplecticForm, dims::Dims{2}) where {T}
    A = zeros(T, dims)
    v = T(J.λ)
    @inbounds for i in 1:size(A, 1)
        if i % 2 == 0
            A[i, i - 1] = -v
        elseif i + 1 ≤ size(A, 2)
            A[i, i + 1] = v
        end
    end
    return A
end
Matrix{T}(J::SymplecticForm, m::Integer, n::Integer) where {T} = Matrix{T}(J, Dims((m, n)))
Matrix(J::SymplecticForm, m::Integer, n::Integer) = Matrix(J, Dims((m, n)))
Array{T}(J::SymplecticForm, dims::Dims{2}) where {T} = Matrix{T}(J, dims)
Array{T}(J::SymplecticForm, m::Integer, n::Integer) where {T} = Matrix{T}(J, m, n)
Array(J::SymplecticForm, m::Integer, n::Integer) = Matrix(J, m, n)
Array(J::SymplecticForm, dims::Dims{2}) = Matrix(J, dims)

Diagonal{T}(J::SymplecticForm, m::Integer) where {T} = Diagonal{T}(fill(zero(T), m))
Diagonal(J::SymplecticForm, m::Integer) = Diagonal{eltype(J)}(J, m)

function dot(A::AbstractMatrix, J::SymplecticForm)
    checksquare(A)
    rows, cols = size(A)
    s = zero(eltype(A))
    @inbounds for i in 1:min(rows, cols)
        if i % 2 == 0
            s += - A[i, i-1]
        elseif i + 1 ≤ cols
            s += A[i, i+1]
        end
    end
    return s * J.λ
end
dot(J::SymplecticForm, A::AbstractMatrix) = dot(A, J)
function dot(x::AbstractVector, J::SymplecticForm, y::AbstractVector) 
    length(x) ≠ length(y) || return dot(x, y)
    s = zero(eltype(x))
    for i in 1:length(x)
        if i % 2 == 0
            @inbounds s += x[i - 1] * J.λ * y[i]
        else
            @inbounds s += -x[i + 1] * J.λ * y[i]
        end
    end
end

muladd(U::UniformScaling, J::SymplecticForm, Z::SymplecticForm) = SymplecticForm(U.λ * J.λ + Z.λ)
muladd(J::SymplecticForm, U::UniformScaling, Z::SymplecticForm) = SymplecticForm(J.λ * U.λ + Z.λ)
muladd(J::SymplecticForm, Z::SymplecticForm, U::UniformScaling) = UniformScaling(- J.λ * Z.λ + U.λ)

_diag_or_value(J::SymplecticForm) = zero(J.λ)

