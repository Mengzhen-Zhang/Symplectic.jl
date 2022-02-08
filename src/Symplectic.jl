module Symplectic

import LinearAlgebra
include("DataTypes/SymplecticForm.jl")
include("DataTypes/CoupledResonators.jl")

export SymplecticForm, 
       Ω,
       CoupledResonators,
       addActive,
       addPassive,
       addGammaEx,
       addGammaIn,
       scatteringMatrix,
       nonSymplecticity, 
       checkSymplectic, 
       modes,
       getQ, 
       getP, 
       ⊗, 
       dsum, 
       ⊕, 
       sct, 
       isct, 
       symplecticCayleyTransform,
       inverseSymplecticCayleyTransform, 
       lstf, 
       localSymplecticTransformation,
       ibs,
       interferenceBasedSequence, 
       decoupleSequence, 
       teleportation, 
       genericity, 
       generify,
       phaseShifting, 
       randomPhaseShifiting, 
       localSymplecticMatrix, 
       randomLocalSymplecticMatrix,
       beamSplitter, 
       circulator, 
       toQPQPBasis, 
       toQQPPBasis, 
       pif,
       preIwasawaFactorization, 
       spectrum

# Utilities
using LinearAlgebra: 
       checksquare, 
       opnorm, 
       eigen, 
       isposdef, 
       issymmetric, 
       norm,
       dot

↑(m, i) = 2*m - (i % 2)
↓(m, i) = 2*m -1 + (i  % 2)

function Lkmij(u, v, m, i, j)
    if norm(v[2*m-1:2*m]) ≈ 0
        return i==j ? 1 : 0
    else
        return (-1)^(j+1)*v[↑(m, i)]*u[↓(m, j)]/( v[2*m]^2 + v[2*m-1]^2 ) + (-1)^i*v[↓(m, i)]*u[↑(m, j)]/( u[2*m]^2 + u[2*m-1]^2 )
    end
end

function getSequence(Ss, modeToDecouple::Int; lambda::Real = 1)
       if length(Ss) == 4
           return decoupleSequence(Ss...; m = modeToDecouple, lambda = lambda)
       else
           step = length(Ss) ÷ 4
           Rs = [ getSequence(Ss[i:i+step-1], modeToDecouple-1 ) for i in 1:step:length(Ss) ]
   
           Ls = decoupleSequence([ *(RLst...) for RLst in Rs]...; m = modeToDecouple, lambda = lambda)
           return vcat([ (i%2==0) ? [ Ls[i] ] : Rs[(i+1) ÷ 2] for i in 1:7 ]...)
       end
end
   
# calculate the 'nonSymplecticity".
nonSymplecticity(J::Union{SymplecticForm, UniformScaling}) = 0
nonSymplecticity(A::AbstractMatrix) = opnorm(transpose(A)*Ω*A - Ω) / opnorm(Ω)

# check if a matrix is symplectic, returning its non-symplecticity
checkSymplectic(J::Union{SymplecticForm, UniformScaling}) = 0
function checkSymplectic(A::AbstractMatrix)
       eltype(A) <: Real || throw(ArgumentError("matrix is not real-valued"))
       checksquare(A) % 2 == 0 || throw(DimensionMismatch("matrix is not even-dimensional"))
       return nonSymplecticity(A)
end
function checkSymplectic(A...)
       nonSymplecticities = []
       for a in A
              if !( a isa Union{SymplecticForm, UniformScaling} )
                     eltype(a) <: Real || throw(ArgumentError("matrix is not real-valued"))
                     checksquare(a) % 2 == 0 || throw(DimensionMismatch("matrix is not even-dimensional"))
              end
              push!(nonSymplecticities, nonSymplecticity(a))
       end
       return nonSymplecticities
end

genericity(S::AbstractMatrix) = (x -> max(x...) / min(x...))(abs.(S))

phaseShifting(θ::Real) = [cos(θ) sin(θ); -sin(θ) cos(θ)]
phaseShifting(θs...) = dsum( map(phaseShifting, θs)... )
randomPhaseShifiting(n::Int) = 
       n == 1 ? phaseShifting(2π * rand()) : phaseShifting( (2π*rand(n))...)

localSymplecticMatrix(θ::Real, l::Real, ϕ::Real) =  phaseShifting(θ) * [l 0 ; 0  1/l ] * phaseShifting(ϕ)
localSymplecticMatrix(θs::Vector, ls::Vector, ϕs::Vector) = 
       dsum( map(args -> localSymplecticMatrix(args...), zip(θs, ls, ϕs))... )
randomLocalSymplecticMatrix(n::Int; bound::Real = 0.2) = 
       n == 1 ? localSymplecticMatrix(2π * rand(), bound*( rand() -1/2) + 1.0, 2π * rand()) :
              localSymplecticMatrix(2π * rand(n), bound*( rand(n) .-1/2) .+ 1.0, 2π * rand(n))

# Generte a Two-Mode BeamSplitter
beamSplitter(angle::Real) = phaseShifting(angle) ⊗ I(2)
# Generate a BeamSplitter between two modes in a multi-mode System
function beamSplitter(angle::Real, mode1, mode2, modes)
    mode1, mode2 = sort([mode1, mode2])
    lst = [ [[2 * (i + 2) - 1, 2 * (i + 2)] for i in 1:mode1-1]... ; [1, 2] ; [[2 * (i + mode1 + 1) - 1, 2 * (i + mode1 + 1)] for i in 1: mode2 - mode1 - 1]...; [3, 4] ; [[2 * (i + mode2 ) - 1, 2 * (i + mode2 )]  for i in 1:modes-mode2]...  ];
    return dsum( beamSplitter(angle), I(2*modes - 4) )[lst, lst]
end

circulator(perm...) = circulator(perm)
circulator(perm::Vector) = 
       hvcat(length(perm), [i == perm[j] ? I(2) : zeros(2, 2)   for i in 1:length(perm), j in 1:length(perm) ]...)

function generify(S::AbstractMatrix; randomize::Function = (n -> I(2 * n)), l::Int = 1)
       n = modes(S)
       seq = vcat( [[randomize(n), S] for i in 1:l ]... , [randomize(n)] )
       g = genericity(*(seq...))
       for i in 1:100
              seqNew = vcat( [[randomize(n), S] for i in 1:l ]... , [randomize(n)] )
              gNew = genericity(*(seqNew...))
              if gNew < g
                     seq = seqNew
                     g = gNew
              end
       end
       return seq
end


modes(A::AbstractMatrix) = ( checkSymplectic(A); size(A, 1) ÷ 2 )

# Get quadtrature (column vector) from a symplectic matrix
getQ(J::Union{SymplecticForm, UniformScaling}, n::Int) = J[:, 2*n-1]
getQ(A::AbstractMatrix, n::Int) = 
       n ≤ modes(A) ? A[:, 2*n-1] : throw(ArgumentError("index is out of bound"))
getQ(J::Union{SymplecticForm, UniformScaling}, n::AbstractRange{<:Integer}) = [J[:, 2*i-1] for i in n]
function getQ(A::AbstractMatrix{T}, n::AbstractRange{<:Integer}) where T
       length(n) ≤ modes(A) || throw(ArgumentError("index is out of bound"))
       return [A[:, 2*i-1] for i in n]
end

getP(J::Union{SymplecticForm, UniformScaling}, n::Int) = J[:, 2*n]
getP(A::AbstractMatrix, n::Int) = 
       n ≤ modes(A) ? A[:, 2*n] : throw(ArgumentError("index is out of bound"))
getP(J::Union{SymplecticForm, UniformScaling}, n::AbstractRange{<:Integer}) = [J[:, 2*i] for i in n]
function getP(A::AbstractMatrix{T}, n::AbstractRange{<:Integer}) where T
       length(n) ≤ modes(A) || throw(ArgumentError("index is out of bound"))
       return [A[:, 2*i] for i in n]
end

const ⊗ = Base.kron

dsum(M1::AbstractMatrix, M2::AbstractMatrix) = cat(M1, M2; dims=(1,2))
dsum(M1::AbstractMatrix, M2::AbstractMatrix, Ms::Vararg{AbstractMatrix}) = dsum(dsum(M1, M2),  mapreduce(identity, dsum, Ms ))
const ⊕ = dsum

# Transforming a symmetric matrix (member of the symplectic algebra) to a symplectic matrix
symplecticCayleyTransform(M::AbstractMatrix) = I - (Ω * M  + I / 2  )^-1
const sct = symplecticCayleyTransform

# Transforming a symplectic matrix to a symmetric matrix
inverseSymplecticCayleyTransform(S::AbstractMatrix) = -Ω*((I - S)^-1 - I/2)
const isct = inverseSymplecticCayleyTransform

# Construct a local symplectic transformation transforming quadrature u to quadrature λΩv (the conjugate of quadrature v)
function localSymplecticTransformation(u::Vector, v::Vector; lambda::Real = 1)
       ( length(u) == length(v) && length(u) % 2 == 0 && length(v) % 2 == 0 ) || throw(ArgumentError("invalid input"))
       n = length(u) ÷ 2
       if n == 1
              return [ Lkmij(u, lambda*v, 1, i, j) for i in 1:2, j in 1:2 ]
       else
              return ⊕([ [ Lkmij(u, lambda*v, m, i, j) for i in 1:2, j in 1:2 ] for m in 1:n ]...)
       end
end
const lstf = localSymplecticTransformation

# Given symplectic matrices S₄, S₃, S₂, and S₁, return a sequence S₄ L₃ S₃ L₂ S₂ L₁ S₁ decoupling the m-th mode, where the local operations are determined up to the scaling factor λ
function decoupleSequence(S4, S3, S2, S1; m::Int=1, lambda::Real = 1)
       L1 = lstf(S1[:,2*m-1], S2[2*m-1,:]; lambda = lambda)
       L3 = lstf(S3[:,2*m-1], S4[2*m-1,:]; lambda = lambda^-1)
       T1 = S2 * L1 * S1
       T2 = S4 * L3 * S3
       L2 = ⊕([ (i==m ? Ω(2) : lstf(T1[2*i-1:2*i,2*m], T2[2*m, 2*i-1:2*i])) for i in 1:(size(T1, 1)÷2)]...)
       return [S4, L3, S3, L2, S2, L1, S1]
end


# Get local Operations L16, L15, ..., L1 for the interfernce-basded sequence LR S L15 S ... S L1 S = ST, with a given target Symplectic Matrix ST
function interferenceBasedSequence(S::AbstractMatrix; T = I(4), lambda::Real = 1)
       return interferenceBasedSequence([S for i in 1:4^modes(T)]; T = T, lambda = lambda)
end
function interferenceBasedSequence(Ss; T = I(4), lambda::Real = 1)    
       n =  modes(Ss[1]) - modes(T)
       nT = modes(T)
       Ss[end] =  Ss[end] * ( inv(T) ⊕ I( 2 * n  ) )
       Ss = getSequence(Ss, nT; lambda = lambda)
       R = ( *(Ss...)[1:2 * nT, 1:2 * nT])^-1  ⊕ I( 2 * n )
       Ss[end] =  Ss[end] * ( T ⊕ I( 2 * n ) )
       return vcat([ R ] , Ss)
end
const ibs = interferenceBasedSequence

function teleportation(S::AbstractMatrix, inModes::Vector, outModes::Vector)
       allModes = [1:(size(S, 1) ÷ 2);]
       In = vcat( [ [2*i - 1, 2 * i] for i in inModes ]... )
       Out = vcat( [ [2*i - 1, 2 * i] for i in outModes ]... )
       ancModes = [i for i in allModes if !(i in inModes)]
       idlModes = [i for i in allModes if !(i in outModes)]
       Usq = 2*ancModes
       Hm = 2*idlModes
       SOutIn = S[ Out, In  ]
       SHmIn = S[ Hm, In ]
       SOutUsq = S[ Out, Usq ]
       SHmUsq = S[ Hm, Usq ]
       return SOutIn - SOutUsq * SHmUsq^-1 * SHmIn
end

function toQQPPBasis(S::AbstractMatrix)
        n = modes(S)
        order = [ [1:2:2*n;]... , [2:2:2*n;]... ]
        return S[order, order]
end

function toQPQPBasis(S::AbstractMatrix)
       n = modes(S)
       order = vcat([ [i, i + n] for i in 1:n]... )
       return S[order, order]
end

function preIwasawaFactorization(S::AbstractMatrix)
       n = modes(S)
       S = toQQPPBasis(S)
       A, B = S[1:n, 1:n], S[1:n, 1+n:2*n]
       C, D = S[1+n:2*n, 1:n], S[1+n:2*n, 1+n:2*n]
       AAtBBt = A*transpose(A) + B*transpose(B)
       P = (C*transpose(A) + D*transpose(B))* AAtBBt^-1
       L = AAtBBt^(1/2)
       X = AAtBBt^(-1/2) * A
       Y = AAtBBt^(-1/2) * B
       return map(toQPQPBasis, [ 
                     [ I(n) zeros(n, n); P I(n) ],
                     [ L zeros(n, n); zeros(n, n) L^-1 ],
                     [ X Y; -Y X]
              ])
end
const pif = preIwasawaFactorization

# Return Symplectic Spectrum Diagonal([λ..., λ...]) and a symplectic matrix S such that Sᵀ M S = toQPQPBasis( Diagonal([λ..., λ...]) )
function spectrum(M::AbstractMatrix)
       (eltype(M)<:Real && isposdef(M)) || throw(ArgumentError("only implemented for real-valued positive-definite matrix"))
       n = size(M, 1) ÷ 2
       K = toQQPPBasis(Ω * M)
       eig = eigen(K; sortby = (λ ->(abs(λ) , imag(λ))))
       spec = [abs(imag( v )) for v in eig.values]
       es = [ (v=eig.vectors[:, 2*i - 1];real(v)) for i in 1:n]
       fs = [ (v=eig.vectors[:, 2*i - 1];-imag(v)) for i in 1:n]
       for i in 1:n
              λ = es[i] ⋅ (toQQPPBasis(Ω(2*n)) * fs[i]) 
              es[i] /= sqrt(λ)
              fs[i] /= sqrt(λ)
       end
       S = hcat( es..., fs... )
       return [spec, toQPQPBasis(S)]
end


end