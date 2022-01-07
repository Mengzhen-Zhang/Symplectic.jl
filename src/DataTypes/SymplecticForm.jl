struct SymplecticForm{T<:Number} 
    λ::T
end

const Ω = SymplecticForm(1)

function (Ω::UniformScaling)(n::Integer) 
    B = zeros(Integer, n, n)
    for i in 1:n
        if i % 2 == 0
            @inbounds B[i, i-1] = -1
        else
            @inbounds B[i, i+1] = 1
        end
    end
    return B
end

Base.eltype(::Type{SymplecticForm{T}}) where {T} = T

Base.ndims(J::SymplecticForm) = 2

Base.has_offset_axes(::SymplecticForm) = false

function Base.getindex(J::SymplecticForm, i::Integer, j::Integer)
    if i % 2 == 0 && j == i - 1 
        return -J.λ 
    elseif i % 2 == 1 && j == i + 1 
        return J.λ 
    else 
        return zero(J.λ)
    end
end
Base.getindex(J::SymplecticForm, n::Integer, m::AbstractRange{<:Integer}) = -getindex(x, m, n)
function Base.getindex(J::SymplecticForm{T}, n::AbstractRange{<:Integer}, m::Integer) where T
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
function Base.getindex(J::SymplecticForm{T}, n::AbstractRange{<:Integer}, m::AbstractRange{<:Integer}) where T
    A = zeros(T, length(n), length(m))
    @inbounds for (j, jj) in enumerate(m), (i, ii) in enumerate(n)
        if ii == jj - 1 && ii % 2 == 0
            v[i] = -J.λ
        elseif ii == jj + 1 && ii % 2 == 1
            v[i] = J.λ
        else
            v[i] = zero(J.λ)
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

function LinearAlgebra.det(J::SymplecticForm{T}) where T
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

@inline LinearAlgebra.mul!(C::AbstractMatrix, A::AbstractMatrix, J::SymplecticForm, alpha::Number, beta::Number) =
    mul!(C, I(size(A, 1)),  A*J, alpha, beta)
@inline LinearAlgebra.mul!(C::AbstractMatrix, J::SymplecticForm, B::AbstractVecOrMat, alpha::Number, beta::Number) = 
    mul!(C, J*B, I(size(B, 2)), alpha, beta)
function LinearAlgebra.mul!(out::AbstractMatrix{T}, a::Number, J::SymplecticForm, α::Number, β::Number) where {T}
    checksquare(out)
    if iszero(β)
        fill!(out, zero(T))
    elseif !isone(β)
        rmul!(out, β)
    end
    s = convert(T, a*J.λ*α)
    if !iszero(s)
        @inbounds for i in diagind(out, -1)
            out[i] += s
        end
        @inbounds for i in diagind(out, 1)
            out[i] -= s
        end
    end
    return out
end
@inline LinearAlgebra.mul!(out::AbstractMatrix{T}, J::SymplecticForm, b::Number, α::Number, β::Number) = 
    mul!(out, J.λ, SymplecticForm(b), α, β)
LinearAlgebra.rmul!(A::AbstractMatrix, J::SymplecticForm) = mul!(A, A, J)
LinearAlgebra.lmul!(J::SymplecticForm, B::AbstractVecOrMat) = mul!(B, J, B)
LinearAlgebra.rdiv!(A::AbstractMatrix, J::SymplecticForm) = mul!(A, A, inv(J))
LinearAlgebra.ldiv!(J::SymplecticForm, B::AbstractVecOrMat) = mul!(B, inv(J), B)
LinearAlgebra.ldiv!(Y::AbstractVecOrMat, J::SymplecticForm, B::AbstractVecOrMat) = (Y .= J \ B)

Broadcast.broadcasted(::typeof(*), x::Number, J::SymplecticForm) = SymplecticForm(x * J.λ)
Broadcast.broadcasted(::typeof(*), J::SymplecticForm, x::Number) = SymplecticForm(J.λ * x)

Broadcast.broadcasted(::typeof(/), J::SymplecticForm, x::Number) = SymplecticForm(J.λ / x)

Broadcast.broadcasted(::typeof(\), x::Number, J::SymplecticForm) = SymplecticForm(x \ J.λ)

function (Base.:^)(J::SymplecticForm, x::Int)
    if x >= 0
        return x % 2 == 0 ? UniformScaling((-J.λ)^(x ÷ 2)) : (-J.λ)^(x ÷ 2) * J
    else
        return inv(J) ^ (- x)
    end
end
Broadcast.broadcasted(::typeof(^), J::UniformScaling, x::Number) = SymplecticForm(J.λ ^ x)

==(J1::SymplecticForm, J2::SymplecticForm) = (J1.λ == J2.λ)
function ==(A::AbstractMatrix, J::SymplecticForm)
    try
        checksquare(A)
    catch e
        return false
    end
    if iszero(J.λ) return iszero(A) end
    return A == J(size(A, 1))
end
==(J::SymplecticForm, A::AbstractMatrix) = A == J

Base.isequal(A::AbstractMatrix, J::SymplecticForm) = false
Base.isequal(J::SymplecticForm, A::AbstractMatrix) = false

function Base.isapprox(J1::SymplecticForm{T}, J2::SymplecticForm(S); 
                atol::Real=0, rtol::Real=Base.rtoldefault(T, S, atol), nans::Bool=false) where {T<:Number, S<:Number}
    return isapprox(J1.λ, J2.λ, rtol=rtol, atol=atol, nans=nans)
end
function Base.isapprox(J::SymplecticForm, A::AbstractMatrix;
                atol::Real=0,
                rtol::Real=Base.rtoldefault(promote_leaf_eltypes(A), eltype(J), atol),
                nans::Bool=false, norm::Function=norm
                )
    n = checksquare(A)
    normJ = norm === opnorm ? abs(J.λ) : 
            norm === LinearAlgebra.norm ? abs(J.λ) * sqrt(n) :
            norm(cat(fill([zero(J.λ) J.λ; -J.λ zero(J.λ)], n)...; dims=(1,2)))
    return norm(A - J) <= max(atol, rtol*max(norm(A), normJ))
end
Base.isapprox(A::AbstractMatrix, J::SymplecticForm; kwargs...) = isapprox(J, A; kwargs...)

function LinearAlgebra.copyto!(A::AbstractMatrix, J::SymplecticForm)
    fill!(A, 0)
    λ = J.λ
    for i = 1:min(size(A, 1), size(A, 2)) - 1
        @inbounds A[i+1, i] = -λ
        @inbounds A[i, i+1] = λ
    end
    return A
end

function LinearAlgebra.cond(J::SymplecticForm{T}) where T
    onereal = inv(one(real(J.λ)))
    return J.λ ≠ zero(T) ? onereal : oftype(onereal, Inf)
end

LinearAlgebra.promote_to_arrays_(n::Int, ::Type{Matrix}, J::SymplecticForm{T}) where{T} = copyto!(Matrix{T}(undef, n, n), J)
LinearAlgebra.promote_to_array_type(J::Tuple{Vararg{AbstractVecOrMat, UniformScaling, SymplecticForm}}) = Matrix

for (f, dim, name) in ((:hcat, 1, "rows"), (:vcat, 2, "cols"))
    @eval begin
        function Base.$f(A::Union{AbstractVecOrMat, UniformScaling, SymplecticForm}...)
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
            n == -1 && throw(ArgumentError($("$f of only Symplectic or UniformScaling objects cannot determine the matrix size")))
            return Base.$f(promote_to_arrays(fill(n, length(A)), 1, promote_to_array_type(A), A...)...)
        end
    end
end

function Base.hvcat(rows::Tuple{Vararg{Int}}, A::Union{AbstractVecOrMat,UniformScaling, SymplecticForm}...)
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

function Base.Matrix{T}(J::SymplecticForm, dims::Dims{2}) where {T}
    A = zeors(T, dims)
    v = T(J.λ)
    for i in diagind(dims..., -1)
        @inbounds A[i] = v
    end
    for i in diagind(dims..., 1)
        @inbounds A[i] = -v
    end
    return A
end
Base.Matrix{T}(J::SymplecticForm, m::Integer, n::Integer) where {T} = Matrix{T}(J, Dims((m, n)))
Base.Matrix(J::SymplecticForm, dims::Dim{2}) = Matrix{eltype(J)}(J, dims)
Base.Matrix{J::SymplecticForm, m::Integer, n::Integer} = Matrix(J, Dims((m, n)))

Base.Array{T}(J::SymplecticForm, dims::Dim{2}) where {T} = Matrix{T}(J, dims)
Base.Array{T}(J::SymplecticForm, m::Integer, n::Integer) where {T} = Matrix{T}(J, m, n)
Base.Array(J::SymplecticForm, m::Integer, n::Integer) = Matrix(J, m, n)
Base.Array(J::SymplecticForm, dims::Dims{2}) = Matrix(J, dims)

LinearAlgebra.Diagonal{T}(J::SymplecticForm, m::Integer) where {T} = Diagonal{T}(fill(zero(T), m))
LinearAlgebra.Diagonal(J::SymplecticForm, m::Integer) = Diagonal{eltype(J)}(J, m)

function LinearAlgebra.dot(A::AbstractMatrix, J::SymplecticForm) 
    rows, cols = size(A)
    s = zero(eltype(A))
    for i in 1:min(rows, cols)
        if i % 2 == 0
            @inbounds s += A[i, i-1]
        elseif i + 1 ≤ cols
            @inbounds s += -A[i, i+1]
        end
    end
    return s * J.λ
end
LinearAlgebra.dot(J::SymplecticForm, A::AbstractMatrix) = dot(A, J)
function LinearAlgebra.dot(x::AbstractVector, J::SymplecticForm, y::AbstractVector) 
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

Base.muladd(J::SymplecticForm, U::UniformScaling, Z::SymplecticForm) = SymplecticForm(J.λ * U.λ + Z.λ)

