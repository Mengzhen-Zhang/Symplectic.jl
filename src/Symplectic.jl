using LinearAlgebra

module Symplectic

using LinearAlgebra

export isSquare, sympCayleyTransform, inverseSympCayleyTransform, sympRound, Symp
export nModes, isGeneric, Omega, Ω, Id, ⊗, dsum, ⊕, randomSymmetric, randomSymp, randomGenericSymp
export monoSymp, localSympFromQuad, getDecoupleSequence, getInterferenceBasedSequence, sympToGraph
export getColors, contract, beamSplitter, circulator

# Used to measure how close a quanity is to zero
tolerance = 10^-6

#=
    Functions not involving Symp
=#
# Check if a Matrix is square
isSquare = LinearAlgebra.checksquare
# Symplectic Cayley Transform: Output a symplectic matrix when the input is a symmetric matrix
function sympCayleyTransform(M::AbstractMatrix)::AbstractMatrix
    n = isSquare(M) ÷ 2
    Ω = cat(fill([0 1; -1 0], n)...; dims=(1,2))
    return I - (Ω * M  + I / 2  )^-1
end
# Invese Symplectic Cayley Transform: Output a symmetric matrix when the input is a symplectic matrix
function inverseSympCayleyTransform(S::AbstractMatrix)::AbstractMatrix
    n = isSquare(S) ÷ 2
    Ω = cat(fill([0 1; -1 0], n)...; dims=(1,2))
    return Symmetric(-Ω*((I - S)^-1 - I/2))
end
# Rounding a Matrix using Cayley Transform
function sympRound(S::AbstractMatrix)::AbstractMatrix
    n = isSquare(S) ÷ 2
    Ω = cat(fill([0 1; -1 0], n)...; dims=(1,2))
    if det(S - I) < tolerance
        return -Ω * ( (sympCayleyTransform ∘ inverseSympCayleyTransform)(Ω * S) )
    else
        return (sympCayleyTransform ∘ inverseSympCayleyTransform)( S )
    end
end


# Data Type For Symplectic Matrices
struct Symp{T}
    S:: AbstractMatrix{T}
    function Symp{T}(S::AbstractMatrix{T}) where T
        if !(T<:Real) return error("not real") end
        n = isSquare(S) ÷ 2
        Ω = cat(fill([0 1; -1 0], n)...; dims=(1,2))
        if (norm(S*Ω*S' - Ω) >= tolerance) return error("not symplectic") end
        return new{T}(S)
    end
end

#=
   Methods
=#

#=
    Property
=#
# Number of Modes
nModes(S::Symp)::Int = size(S.S)[1] ÷ 2
# Symmetric?
LinearAlgebra.issymmetric(S::Symp)::Bool = issymmetric(S.S)
# Positive-Definite?
LinearAlgebra.isposdef(S::Symp)::Bool = isposdef(S.S)
# Diagonal?
LinearAlgebra.isdiag(S::Symp)::Bool = isdiag(S.S)

# Generic?
isGeneric(m::AbstractMatrix)::Bool = all(x -> abs(x) > tolerance, m )
isGeneric(S::Symp)::Bool = isGeneric(S.S)

#=
    Constructor
=#    
# From Matrix
Symp(S) = Symp{typeof(S).parameters[1]}(S)
# From Symplectic Matrix
Symp(S::Symp) = S
# Rounding a Symplectic Matrix using Cayley Transform
sympRound(S::Symp)::Symp = sympRound(S.S)
# From Number of Modes to Symplectic Form
Omega(n::Int)::Symp = cat(fill([0 1; -1 0], n)...; dims=(1,2)) 
Ω = Omega
# From Number of Modes to Identity Symplectic Matrix
Id(n::Int)::Symp = Matrix(1I, 2*n, 2*n)

#= 
    Conversion
=#
Base.convert(::Type{Symp}, S) = Symp(S)
Base.convert(::Type{Symp}, S::Symp) = S
Base.convert(::Type{AbstractMatrix}, S::Symp) = S.S
Base.Matrix(S::Symp) = S.S


#= 
    IO
=#
Base.print(S::Symp) = print(S.S)
Base.println(S::Symp) = println(S.S)
Base.show(S::Symp) = show(S.S)
Base.show(io::IO, ::MIME"text/plain", S::Symp) = show(io, MIME"text/plain"(), S.S)

#=
    Unary Operators
=#
# Transpose
Base.transpose(S::Symp)::Symp = Base.transpose( S.S )
Base.adjoint(S::Symp)::Symp = Base.adjoint( S.S )
# Unary Minus and Unary Plus
Base.:-(S::Symp)::Symp = Base.:-( S.S )
Base.:+(S::Symp)::Symp = Base.:+( S.S )
# Inverse
Base.inv(S::Symp)::Symp = - Ω(nModes(S)) * S' * Ω(nModes(S)) 
# Norm
LinearAlgebra.norm(S::Symp)::Number = norm(S.S)

#=
    Binary Operators
=#

# Multiplication
Base.:*(S::Symp, M::AbstractMatrix)::AbstractMatrix = *( S.S,  M )
Base.:*(M::AbstractMatrix, S::Symp)::AbstractMatrix = *( M,  S.S )
Base.:*(S1::Symp, S2::Symp)::Symp = *( S1.S,  S2.S )
Base.:*(p::Number, S::Symp)::AbstractMatrix = *(p, S.S)
Base.:*(S::Symp, p::Number)::AbstractMatrix = *(S.S, p)
Base.:*(S1::Symp, S2::Symp, Ss::Vararg{Symp})::Symp = *(S1.S, S2.S,  [s.S for s in Ss]... )
# Addition
Base.:+(S::Symp, M::AbstractMatrix)::AbstractMatrix =  +( S.S,  M )
Base.:+(M::AbstractMatrix, S::Symp)::AbstractMatrix =  +( M,  S.S )
Base.:+(S1::Symp, S2::Symp)::AbstractMatrix = +( S1.S,  S2.S )
# Minus
Base.:-(S::Symp, M::AbstractMatrix)::AbstractMatrix = -( S.S,  M )
Base.:-(M::AbstractMatrix, S::Symp)::AbstractMatrix = -( M,  S.S )
Base.:-(S1::Symp, S2::Symp)::AbstractMatrix = -( S1.S,  S2.S )
# Tensor Product
Base.kron(S::Symp, M::AbstractMatrix)::AbstractMatrix = kron( S.S,  M )
Base.kron(M::AbstractMatrix, S::Symp)::AbstractMatrix = kron( M,  S.S )
Base.kron(S1::Symp, S2::Symp) = kron(S1.S, S2.S )
⊗ = Base.kron
# Direct Sum
dsum(M1::AbstractMatrix, M2::AbstractMatrix)::AbstractMatrix = cat(M1, M2; dims=(1,2))
dsum(M1::AbstractMatrix, M2::AbstractMatrix, Ms::Vararg{AbstractMatrix})::AbstractMatrix = dsum(dsum(M1, M2),  mapreduce(identity, dsum, Ms ))
dsum(S::Symp, M::AbstractMatrix)::AbstractMatrix = dsum( S.S,  M )
dsum(M::AbstractMatrix, S::Symp)::AbstractMatrix = dsum( M,  S.S )
dsum(S1::Symp, S2::Symp)::Symp = dsum( S1.S,  S2.S )
dsum(S1::Symp, S2::Symp, Ss::Vararg{Symp})::Symp = dsum(S1.S, S2.S,  [s.S for s in Ss]... )
⊕ = dsum
# Matrix Power
Base.:^(S::Symp, p::Int)::Symp = p < 0 ? inv(S).S ^ (-p)  : S.S^p
Base.:^(S::Symp, p::Number)::AbstractMatrix = S.S^p
# Exponetiation
Base.:^(p::Number, S::Symp)::AbstractMatrix = p^S.S
Base.exp(S::Symp)::AbstractMatrix = exp(S.S)
# Division
Base.:/(S::Symp, p::Number)::AbstractMatrix = S.S / p

# Functions

# From Number of Modes to Random Symmetric Matrix
function randomSymmetric(n::Int, range::Real = 3)::Matrix{BigFloat}
    if range <= 0
        return error("range must be positive")
    end
    return Symmetric( range*(rand(2*n, 2*n) - I ) )
end
# From Number of Modes to Random Symplectic Matrix
function randomSymp(n::Int, range::Real = 3)::Symp
    S1 = Id(n) - (Ω(n) * randomSymmetric(n, range) + Id(n) / 2  )^-1
    S2 = Id(n) - (Ω(n) * randomSymmetric(n, range) + Id(n) / 2  )^-1
    return  S1 * S2 
end
function randomGenericSymp(n::Int, range::Real = 3)::Symp
        m = randomSymp(n, range)
        return isGeneric(m) ? m : randomGenericSymp(n, range)
end
# Single Mode Gaussian Operation
monoSymp(θ::Real)::Symp =  exp(θ*Ω(1)) 
# From (source quadrautre, target quadratre) to Local Symplectic Matrix transforming source to Ω target
↑(m, i) = 2*m - (i % 2)
↓(m, i) = 2*m -1 + (i  % 2)
getPattern(v::Vector)=(n->n>tolerance).( (diag∘sqrt∘Diagonal∘(x->x'*x)∘(v->reshape(v, 2, :)))(v) )
function Lkmij(u, v, m, i, j)
    if norm(v[2*m-1:2*m])<tolerance
        return i==j ? 1 : 0
    else
        return (-1)^(j+1)*v[↑(m, i)]*u[↓(m, j)]/( v[2*m]^2 + v[2*m-1]^2 ) + (-1)^i*v[↓(m, i)]*u[↑(m, j)]/( u[2*m]^2 + u[2*m-1]^2 )
    end
end
function localSympFromQuad(u::Vector, v::Vector)::Symp
    if length(u) % 2 !=0 || length(u) != length(v) || getPattern(u) != getPattern(v) 
        return error("input not valid") 
    end
    n = length(u) ÷ 2
    return ⊕([ Symp([ Lkmij(u, v, m, i, j) for i in 1:2, j in 1:2 ]) for m in 1:n ]...)
    return mapreduce( identity, ⊕, 
        [ Symp([ Lkmij(u, v, m, i, j) for i in 1:2, j in 1:2 ]) for m in 1:n ])
end
# Get Local Operations L3 L2 L1 for the Decoupling Sequence S4 L3 S3 L2 S2 L1 S1
function getDecoupleSequence(S4::Symp, S3::Symp, S2::Symp, S1::Symp, m::Int=1)::Vector{Symp}
    # if !(all(x->isGeneric(x), [S1, S2, S3, S4]) ) return error("not generic") end
    L1 = localSympFromQuad(S1.S[:,2*m-1], S2.S[2*m-1,:])
    L3 = localSympFromQuad(S3.S[:,2*m-1], S4.S[2*m-1,:])
    T1 = S2 * L1 * S1
    T2 = S4 * L3 * S3
    u = T1.S[:,2*m]
    v = T2.S[2*m,:] + sum( (k->T1.S[k,2*m-1]*T1.S[k,2*m]-T2.S[2*m-1,k]*T2.S[2*m,k]).(1:length(u))  ) * T2.S[2*m-1, :]
    L2 = localSympFromQuad(u, v)
    return [S4, L3, S3, L2, S2, L1, S1]
end
# Get local Operations L16, L15, ..., L1 for the interfernce-basded sequence LR S L15 S ... S L1 S, with a given 4×4 target Symplectic Matrix S⊙
function getSequence(Ss::Vector{Symp{T}}, modeToDecouple::Int)::Vector{Symp} where T
    if length(Ss) == 4
        return getDecoupleSequence(Ss..., modeToDecouple)
    else
        step = length(Ss) ÷ 4
        Rs = [ getSequence(Ss[i:i+step-1], modeToDecouple-1 ) for i in 1:step:length(Ss) ]

        Ls = getDecoupleSequence([ *(RLst...) for RLst in Rs]..., modeToDecouple)
        return vcat([ (i%2==0) ? Ls[i] : Rs[(i+1) ÷ 2] for i in 1:7 ]...)
    end
end
function getInterferenceBasedSequence(S::Symp, ST::Symp)::Vector{Symp}
    if !(isGeneric(S))  return error("not generic") end
    Ss = [S for i in 1:4^nModes(ST)]
    Ss[end] =  Ss[end] * (ST ⊕ Id(nModes(S) - nModes(ST)))^-1
    Ss = getSequence(Ss, nModes(ST))
    LR = Symp( ( *(Ss...).S[1:2*nModes(ST), 1:2*nModes(ST)]) )^-1  ⊕ Id(nModes(S) - nModes(ST)) 
    Ss[end] =  Ss[end] * (ST ⊕ Id(nModes(S) - nModes(ST)))
    return vcat([ LR ] , Ss...)
end
# Generate a Graph based on the given Symplectic Matrix
const Graph = AbstractMatrix{Bool}
const Colors = Vector
const Vertex = Int
const Colored = Bool

function sympToGraph( S::Symp )::Graph
    n = nModes(S)
    return [ norm(S.S[2*i-1:2*i, 2*j-1:2*j]) > tolerance for i in 1:n, j in 1:n ]
end

function colorSuccessors(sucs, Cs::Colors)::Colors
    firstColor = Cs[sucs[1]]
    Cs[sucs] .= firstColor
    return Cs
end  

function  getColors( G::Graph, Cs::Colors )::Colors
    colored = false
    n = length(Cs)
    for v in 1:n
        sameColorVs = G[ findall(x->x == Cs[v], Cs ), :]
        sucs = findall(vi->any(sameColorVs[:, vi] ),  [1:n;] )
        nColors = length( Set(Cs[sucs]) )
        if nColors > 1
            Cs = colorSuccessors( sucs, Cs )
            colored = true
        end
    end
    return colored ? getColors(G, Cs) : Cs
end

# Is i a successor of j (both colored)?
function areColoredVeticesConnected(ci, cj, colorSets::Dict, G::Graph)::Bool
    vis = colorSets[ci]
    vjs = colorSets[cj]
    return any(G[vis, vjs])
end

function contract( G::Graph, Cs::Colors )
    allColors = [Set(Cs)...]
    colorSets = Dict(c => findall(x -> x == c, Cs) for c in allColors)
    return [areColoredVeticesConnected(ci, cj, colorSets, G) for  ci in allColors, cj in allColors]
end

# Generte a Two-Mode BeamSplitter
beamSplitter(angle::Real)::Symp = monoSymp(angle) ⊗ Id(1)

# Generate a Circulator-Like Symplectic Matrix according to the given permutation
circulator(perm::Vector)::Symp = hvcat(length(perm), [i == perm[j] ? Id(1).S : zeros(2, 2)   for i in 1:length(perm), j in 1:length(perm) ]...)

end