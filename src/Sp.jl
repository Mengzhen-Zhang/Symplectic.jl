abstract type AbstractSp <: AbstractFunc end

struct Sp <: AbstractSp
    l::Integer
    f::Function
    dim::integer
end

convert(Func, sp::AbstractSp) = Func(get_l(sp), get_f(sp))

convert(Sp, sp::AbstractSp) = Sp(get_l(sp), get_f(sp), get_dim(sp))

get_l(sp::Sp) = sp.l            # Interface
get_f(sp::Sp) = sp.f            # Interface
get_dim(sp::Sp) = sp.dim        # Interface


Ω(sp::AbstractSp) = begin
    d = get_dim(sp)
    Sp(0, Ω(2 * d), d)
end

Base.one(sp::AbstractSp) = begin
    d = get_dim(sp)
    Sp(0, I(2 * d), d)
end

Base.inv(sp::AbstractSp) = - Ω(sp) * transpose(sp) * Ω(sp)

extend(sp::AbstractMatrix, n::Integer) = dsum(sp, Func(0, (xs::Vararg) -> I(2*n)))
extend(n::Integer, sp::AbstractMatrix) = dsum(Func(0, (xs::Vararg) -> I(2*n)), sp)


# Beamsplitter

struct BeamsplitterSp
    l::Integer
    i::Integer
    j::Integer
    phase::Function
end

get_l(sp::BeamsplitterSp) = sp.l
get_i(sp::BeamsplitterSp) = sp.i
get_j(sp::BeamsplitterSp) = sp.j
get_phase(sp::BeamsplitterSp) = sp.phase

get_dim(sp::BeamsplitterSp) = max(get_i(sp),get_j(sp))

function convert(Sp, sp::BeamsplitterSp)
    l = get_l(sp)
    dim = get_dim(sp)
    i, j = get_i(sp), get_j(sp)
    phase = get_phase(sp)
    f(xs::Vararg) = beamSplitter(phase(xs...), i, j, dim)
    return Sp(l, i, j phase)
end

function Base.transpose(sp::BeamsplitterSp)
    l, i, j, phase = get_l(sp), get_i(sp), get_j(sp), get_phase(sp)
    return BeamsplitterSp(l, j, i, phase)
end

Base.inv(sp::BeamsplitterSp) = Base.transpose(sp::BeamsplitterSp)


# Amplifier

struct AmpliferSp
    l::Integer
    i::Integer
    j::Integer
    amplification::Function
end

get_l(sp::AmpliferSp) = sp.l
get_i(sp::AmpliferSp) = sp.i
get_j(sp::AmpliferSp) = sp.j
get_amplification(sp::AmpliferSp) = sp.amplification

get_dim(sp::AmpliferSp) = max(get_i(sp),get_j(sp))

function convert(Sp, sp::AmpliferSp)
    l = get_l(sp)
    dim = get_dim(sp)
    i, j = get_i(sp), get_j(sp)
    amplification = get_amplification(sp)
    f(xs::Vararg) = Amplifer(amplification(xs...), i, j, dim)
    return Sp(l, i, j amplification)
end

function Base.inv(sp::AmpliferSp)
    l, i, j, amplification = get_l(sp), get_i(sp), get_j(sp), get_amplification(sp)
    return AmpliferSp(l, i, j, amplification)
end

Base.transpose(sp::AmpliferSp) = sp






