struct MatrixSequence{T} # <:  Vector{AbstractMatrix{T}}
    sequence:: Array{Array{T,2},1}
    interspersion:: Array{Array{T,2},1}
    size::Int
    function MatrixSequence{T}(sequence::Array{Array{T,2},1}, interspersion::Array{Array{T,2},1}) where T
        if (length(sequence) < 1)
            return error("not a sequence")
        end
        if (length(interspersion) != length(sequence) + 1 )
            return error("not an interspersion")
        end
        if size(sequence[1])[1] != size(interspersion[1])[1]
            return error("inconsitent sizes")
        end
        for matrix in sequence
            if !(LinearAlgebraUtilities.isSquare(matrix))
                return error("not square")
            end
            if size(matrix) != size(sequence[1])
                return error("inconsitent sizes")
            end
        end
        for matrix in interspersion
            if !(LinearAlgebraUtilities.isSquare(matrix))
                return error("not square")
            end
            if size(matrix) != size(sequence[1])
                return error("inconsitent sizes")
            end
        end
        
        return new{T}(
            map(AbstractMatrix{T}, sequence),
            map(AbstractMatrix{T}, interspersion),
            size(sequence[1])[1]
        )
    end
end
function MatrixSequence(sequence::Array{Array{T,2},1}, interspersion) where T
    return MatrixSequence{T}(sequence, Array{Array{T,2}, 1}(interspersion))
end
function MatrixSequence(sequence::Array{Array{T,2}, 1}) where T
    return MatrixSequence{T}(sequence, Array{Array{T,2}, 1}([LinearAlgebra.I(size(sequence[1])[1]) for i in 1:size(sequence)[1]+1]))
end
function MatrixSequence(matrixSequence::MatrixSequence)
    matrixSequence
end
function Base.convert(::Type{MatrixSequence}, matrixSequence::MatrixSequence)
    matrixSequence
end
function Base.convert(::Type{AbstractMatrix}, matrixSequence::MatrixSequence)
    Matrix(matrixSequence)
end
function Base.size(matrixSequence::MatrixSequence, dim::Integer)
    size(matrixSequence.sequence, dim)
end
function Base.size(matrixSequence::MatrixSequence)
    size(matrixSequence.sequence)
end
function Base.length(matrixSequence::MatrixSequence)
    return length(matrixSequence.sequence)
end
function Base.Matrix(matrixSequence::MatrixSequence)::AbstractMatrix
    sequence = matrixSequence.sequence
    interspersion =  matrixSequence.interspersion
    return foldl(
        (cur, i)-> cur * sequence[i]*interspersion[i+1],
        1:size(matrixSequence)[1];
        init = interspersion[1])
end
function sliceInner(matrixSequence::MatrixSequence, i::Integer, j::Integer)::MatrixSequence
    sequence, interspersion, Size = matrixSequence.sequence, matrixSequence.interspersion, matrixSequence.size
    return MatrixSequence(sequence[i:j], [LinearAlgebra.I(Size),  interspersion[i+1:j]..., LinearAlgebra.I(Size)]) 
end
function slice(matrixSequence::MatrixSequence, i::Integer, j::Integer)::MatrixSequence
    sequence, interspersion, Size = matrixSequence.sequence, matrixSequence.interspersion, matrixSequence.size
    return MatrixSequence(sequence[i:j], interspersion[i:j+1]) 
end
function Base.:*(MS1::MatrixSequence, MS2::MatrixSequence)
    sequence1, inspersion1 = MS1.sequence, MS1.interspersion
    sequence2, inspersion2 = MS2.sequence, MS2.interspersion
    sequence = [sequence1 ; sequence2]
    inspersion = [inspersion1[1:end - 1] ; [ inspersion1[end] * inspersion2[1] ] ; inspersion2[2:end]]
    return MatrixSequence(sequence, inspersion)
end
function setInterspersion(matrix::AbstractMatrix, matrixSequence::MatrixSequence, i::Integer)::MatrixSequence
    sequence = matrixSequence.sequence
    interspersion = matrixSequence.interspersion
    matrixSequence.interspersion[i] = matrix
    return MatrixSequence(sequence, interspersion)
end
function setSequence(matrix::AbstractMatrix, matrixSequence::MatrixSequence, i::Integer)::MatrixSequence
    sequence = matrixSequence.sequence
    interspersion = matrixSequence.interspersion
    matrixSequence.sequence[i] = matrix
    return MatrixSequence(sequence, interspersion)
end
function replaceSequence(subMatrixSequence::MatrixSequence, matrixSequence::MatrixSequence, i::Integer)
    subSequence, subInterspersion = subMatrixSequence.sequence, subMatrixSequence.interspersion
    sequence, interspersion = matrixSequence.sequence, matrixSequence.interspersion
    return MatrixSequence(
        [sequence[1:i-1] ; subSequence ; sequence[i+1:end]],
        [interspersion[1:i-1] ; [interspersion[i]*subInterspersion[1]] ; subInterspersion[2:end-1] ; [subInterspersion[end]*interspersion[i+1]] ; interspersion[i+2:end]]
    )
end