
# Orthogonal Symplectic Operations

abstract type AbstractOSp <: AbstractSp end 

struct OSp <: AbstractOSp
    X::AbstractMatrix
    Y::AbstractMatrix
end

convert(Sp, sp::AbstractOSp) = Sp(Matrix(sp))

get_X(sp::AbstractOSp) = get_X(sp::OSp) # Interface
get_Y(sp::AbstractOSp) = get_X(sp::OSp) # Interface

function Matrix(sp::OSp)
    X, Y = get_X(sp), get_Y(sp)
    return  qpqp([X Y; -Y X])
end

dim(sp::AbstractOSp) = size(get_X(sp), 1)

function Base.transpose(sp::AbstractOSp)
    X, Y = get_X(sp), get_Y(sp)
    transposed_X = transpose(X)
    transposed_Y = transpose(Y)
    minus_transposed_Y = [variable_symbol ? e == variable_symbol : - e for e in transposed_Y]
    return OSp(transposed_X, minus_transposed_Y)
end

Base.inv(sp::AbstractOSp) = transpose(sp)


# Beamsplitter

struct BeamsplitterSp <: AbstractOSp
    phase
    i::Integer
    j::Integer
end

get_i(sp::BeamsplitterSp) = sp.i
get_j(sp::BeamsplitterSp) = sp.j
get_phase(sp::BeamsplitterSp) = sp.phase

dim(sp::BeamsplitterSp) = max(get_i(sp), get_j(sp))

convert(OSp, sp::BeamsplitterSp) = OSp(get_X(sp), get_Y(sp))


function get_X(sp::BeamsplitterSp)
    i, j = get_i(sp), get_j(sp)
    c, d = min(i, j), max(i, j)

    if phase == variable_symbol
        
end


