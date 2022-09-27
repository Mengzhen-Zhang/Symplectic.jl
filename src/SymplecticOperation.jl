export toSympType, SymplecticOperation
export SympType, BeamSplitter, PhaseShifting, Custom, Undef

const Params = AbstractArray{Any}

@enum SympType PhaseShifting BeamSplitter Custom Undef

function toSympType(name::AbstractString)
    if name == "beam_splitter"
        return BeamSplitter
    elseif name == "phase_shifting"
        return PhaseShifting
    elseif name == "custom"
        return Custom
    else
        return Undef
    end
end

function funcBuild(type::SympType, params::Params)::Func
    l = length([0 for p in params if p == "_"])

    func = if type == PhaseShifting
        phaseShifting
    elseif type == BeamSplitter
        beamSplitter
    elseif type == Custom
        customReshape
    else
        throw(ArgumentError("Operation undefined"))
    end

    f = (xs::Vararg) -> begin
        i = 1
        ys = []
        for p in params
            if p == "_"
                push!(ys, xs[i])
                i += 1
            else
                push!(ys, p)
            end
        end
        func(ys...)
    end
    return Func(l, f)
end


struct SymplecticOperation
    n::Int64
    Op::Func
    function SymplecticOperation(n::Int64, Op::Func)
        return new(n, Op)
    end
end

function (s::SymplecticOperation)(xs::Vararg)
    return s.Op(xs...)
end

function SymplecticOperation(n::Int64, type::SympType, params::Params)::SymplecticOperation
    op = funcBuild(type, params)
    return SymplecticOperation(n, op)
end

function SymplecticOperation(m::AbstractMatrix)
    n = size(m, 1) ÷ 2
    l = 0
    f = (xs::Vararg) -> m
    return SymplecticOperation(n, Func(l, f))
end

function dsum(s1::SymplecticOperation, s2::SymplecticOperation)::SymplecticOperation
    n = s1.n + s2.n
    l1, l2 = s1.Op.l, s2.Op.l
    l = l1 + l2
    op = (xs::Vararg) -> begin
        dsum(s1(xs[1:l1]...), s2(xs[l1+1:l]...))
    end
    return SymplecticOperation(n, Func(l, op))
end

function dsum(s1::SymplecticOperation, n2::Int64)::SymplecticOperation
    n = s1.n + n2
    l = s1.Op.l
    op = (xs::Vararg) -> begin
        dsum(s1(xs...), I(2 * n2))
    end
    return SymplecticOperation(n, Func(l, op))
end

function dsum(n1::Int64, s2::SymplecticOperation)::SymplecticOperation
    n = n1 + s2.n
    l = s2.Op.l
    op = (xs::Vararg) -> begin
        dsum(I(2 * n1), s2(xs...))
    end
    return SymplecticOperation(n, Func(l, op))
end

function Base.:*(s1::SymplecticOperation, s2::SymplecticOperation)::SymplecticOperation
    if s1.n < s2.n
        s1 = dsum(s1, s2.n - s1.n)
    end
    if s1.n > s2.n
        s2 = dsum(s2, s1.n - s2.n)
    end

    n = s1.n
    op = s1.Op * s2.Op

    return SymplecticOperation(n, op)
end

function Base.inv(s::SymplecticOperation)
    n = s.n
    l = s.Op.l
    op = (xs::Vararg) -> begin
        -Ω(2 * n) * transpose(s(xs...)) * Ω(2 * n)
    end
    return SymplecticOperation(n, Func(l, op))
end

function teleportation(s::SymplecticOperation, inModes, outModes)::SymplecticOperation
    l = s.Op.l
    n = length(inModes)
    if n == 0
        n = s.n
    end
    op = (xs::Vararg) -> begin
        teleportation(s(xs...), inModes, outModes)
    end
    return SymplecticOperation(n, Func(l, op))
end

function nonSymplecticity(s::SymplecticOperation)::Func
    l = s.Op.l
    f = (xs::Vararg) -> begin
        nonSymplecticity(s(xs...))
    end
    return Func(l, f)
end
