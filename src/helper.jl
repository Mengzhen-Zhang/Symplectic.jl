using LinearAlgebra: I, norm, eigen, checksquare, UniformScaling

# TODO: deprecate
"""
    checkSymplectic(A...)

Check if a sequence of matrices are all symplectic. If a matrix is symplectic,
returning its non-symplecticity
"""
checkSymplectic(J::Union{SymplecticForm,UniformScaling}) = 0
function checkSymplectic(A::AbstractMatrix)
       eltype(A) <: Real || throw(ArgumentError("matrix is not real-valued"))
       checksquare(A) % 2 == 0 || throw(DimensionMismatch("matrix is not even-dimensional"))
       return nonSymplecticity(A)
end
function checkSymplectic(A...)
       nonSymplecticities = []
       for a in A
              if !(a isa Union{SymplecticForm,UniformScaling})
                     eltype(a) <: Real || throw(ArgumentError("matrix is not real-valued"))
                     checksquare(a) % 2 == 0 || throw(DimensionMismatch("matrix is not even-dimensional"))
              end
              push!(nonSymplecticities, nonSymplecticity(a))
       end
       return nonSymplecticities
end

# TODO: deprecate
genericity(S::AbstractMatrix) = (x -> max(x...) / min(x...))(abs.(S))

#TODO: deprecate
function generify(S::AbstractMatrix; randomize::Function=(n -> I(2 * n)), l::Int=1)
       n = modes(S)
       seq = vcat([[randomize(n), S] for i in 1:l]..., [randomize(n)])
       g = genericity(*(seq...))
       for i in 1:100
              seqNew = vcat([[randomize(n), S] for i in 1:l]..., [randomize(n)])
              gNew = genericity(*(seqNew...))
              if gNew < g
                     seq = seqNew
                     g = gNew
              end
       end
       return seq
end

# TODO: deprecate
modes(A::AbstractMatrix) = (checkSymplectic(A); size(A, 1) ÷ 2)

# TODO: deprecate
# Get quadtrature (column vector) from a symplectic matrix
getQ(J::Union{SymplecticForm,UniformScaling}, n::Int) = J[:, 2*n-1]
getQ(A::AbstractMatrix, n::Int) =
       n ≤ modes(A) ? A[:, 2*n-1] : throw(ArgumentError("index is out of bound"))
getQ(J::Union{SymplecticForm,UniformScaling}, n::AbstractRange{<:Integer}) = [J[:, 2*i-1] for i in n]
function getQ(A::AbstractMatrix{T}, n::AbstractRange{<:Integer}) where {T}
       length(n) ≤ modes(A) || throw(ArgumentError("index is out of bound"))
       return [A[:, 2*i-1] for i in n]
end

# TODO: deprecate
getP(J::Union{SymplecticForm,UniformScaling}, n::Int) = J[:, 2*n]
getP(A::AbstractMatrix, n::Int) =
       n ≤ modes(A) ? A[:, 2*n] : throw(ArgumentError("index is out of bound"))
getP(J::Union{SymplecticForm,UniformScaling}, n::AbstractRange{<:Integer}) = [J[:, 2*i] for i in n]
function getP(A::AbstractMatrix{T}, n::AbstractRange{<:Integer}) where {T}
       length(n) ≤ modes(A) || throw(ArgumentError("index is out of bound"))
       return [A[:, 2*i] for i in n]
end

# TODO: deprecate
const ⊗ = Base.kron

# TODO: deprecate
"""
    customReshape(xs::Vararg)

Reshape the sequence xs[1:end-1] into a matrix of shape xs[end] × xs[end].
"""
function customReshape(xs::Vararg)
       return reshape([xs[1:end-1]...], (xs[end], xs[end]))
end
