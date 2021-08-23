export Symp, isSquare, isGeneric, nModes, ⊗, dsum, ⊕

const tolerance = 10^-6

struct Symp{T} <: AbstractMatrix{T}
    S:: AbstractMatrix{T}
    function Symp{T}(S::AbstractMatrix{T}) where T
        if !(T<:Real) return error("not real") end
        n = LinearAlgebra.checksquare(S) ÷ 2
        Ω = cat(fill([0 1; -1 0], n)...; dims=(1,2))
        if (norm(S*Ω*S' - Ω) >= tolerance) return error("not symplectic") end
        return new{T}(S)
    end
end

# Define AbstractMatrix Interface `size`
Base.size(S::Symp, dim::Integer)::Integer = size(S.S, dim)
Base.size(S::Symp) = size(S.S)
# Define AbstractMatrix Interface `getindex`
function Base.getindex(S::Symp{T}, i::Integer, j::Integer)::T where {T} 
    return S.S[i, j]
end
# Type Conversions
Symp(S::AbstractMatrix) = Symp{typeof(S).parameters[1]}(S)
Symp(S::Symp) = S
Base.Matrix(S::Symp)::AbstractMatrix = S.S
Base.convert(::Type{Symp}, S::AbstractMatrix) = Symp(S)
Base.convert(::Type{Symp}, S::Symp) = S
Base.convert(::Type{AbstractMatrix}, S::Symp) = Matrix(S)


# Properties
isSquare = LinearAlgebra.checksquare
LinearAlgebra.issymmetric(S::Symp)::Bool = issymmetric(S.S)
LinearAlgebra.isposdef(S::Symp)::Bool = isposdef(S.S)
LinearAlgebra.isdiag(S::Symp)::Bool = isdiag(S.S)
isGeneric(m::AbstractMatrix)::Bool = all(x -> abs(x) > tolerance, m )
isGeneric(S::Symp)::Bool = isGeneric(S.S)
nModes(S::Symp)::Int = size(S)[1] ÷ 2       # Number of Modes

# Operators
Base.transpose(S::Symp)::Symp = Base.transpose( S.S )
Base.adjoint(S::Symp)::Symp = Base.adjoint( S.S )
# Unary Minus and Unary Plus
Base.:-(S::Symp)::Symp = Base.:-( S.S )
Base.:+(S::Symp)::Symp = Base.:+( S.S )
# Inverse
Base.inv(S::Symp)::Symp = - Ω(nModes(S)) * S' * Ω(nModes(S)) 
Base.:*(S1::Symp, S2::Symp, Ss::Vararg{Symp})::Symp = *(S1.S, S2.S,  [s.S for s in Ss]... )
⊗ = Base.kron
# Direct Sum
dsum(M1::AbstractMatrix, M2::AbstractMatrix)::AbstractMatrix = cat(M1, M2; dims=(1,2))
dsum(M1::AbstractMatrix, M2::AbstractMatrix, Ms::Vararg{AbstractMatrix})::AbstractMatrix = dsum(dsum(M1, M2),  mapreduce(identity, dsum, Ms ))
dsum(S1::Symp, S2::Symp, Ss::Vararg{Symp})::Symp = dsum(S1.S, S2.S,  [s.S for s in Ss]... )
⊕ = dsum
Base.:^(S::Symp, p::Int)::Symp = p < 0 ? inv(S).S ^ (-p)  : S.S^p