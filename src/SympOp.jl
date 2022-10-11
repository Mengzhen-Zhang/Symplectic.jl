# depend on SymplecticForm.jl
# depend on phase_shifting.jl
# depend on beam_splitter.jl
# depend on amplifier.jl

export sympop, SympOp

using LinearAlgebra: I

struct SympOp
    nargs::Int                  # length of arguments of fn
    nmodes::Int                 # number of modes
    fn::Function                # 2nmodes×2nmodes-matrix-valued function taking a vector of length nargs as input
    ensure::Bool                # if true, add a restriction to the cost function
    function SympOp(nargs::Int, nmodes::Int, fn::Function, ensure::Bool)
        if nargs < 0
            return throw(ArgumentError("number of args must be non-negative"))
        end
        if nmodes < 1
            return throw(ArgumentError("number of modes must be more than zero"))
        end
        return new(nargs, nmodes, fn, ensure)
    end
end

"""
    sympop(nargs::Int, nmodes::Int, fn::Function, ensure::Bool=false)

Return a `SympOp` object.

An `SympOp` object is a wrapper of a 2(nmodes)×2(nmodes)-matrix-valued function `fn` which takes
a length `nargs` real vector as the input. Return of `fn` is expected to be symplectic. When it is
not possible to ensure the symplecticity of `fn`'s return, the user can set `ensure=true`, so that
in other functions of this package where symplecticity is required, an automatic check will be
executed.
"""
function sympop(nargs::Int, nmodes::Int, fn::Function; ensure::Bool=false)
    return SympOp(nargs, nmodes, fn, ensure)
end

"""
    sympop(m::AbstractMatrix, ensure::Bool=false)

Construct a sympop from a constant matrix.
"""
function sympop(m::AbstractMatrix, ensure::Bool=false)
    if size(m ,1) != size(m, 2) || size(m, 1) % 2 != 0
        throw(ArgumentError("must be even-dimensional square matrix"))
    end
    return sympop(0, size(m, 1)÷2, (args::AbstractArray)->m; ensure=ensure)
end

function sympop(type::Symbol, params::Vararg)
    # phase-shifting
    if type == :ps              
        if length(params) == 0
            return sympop(1, 1, (args::AbstractArray)->phaseShifting(args[1]))
        elseif length(params) = 1
            return sympop(1, params[1], (args::AbstractArray)->phaseShifting(args...))
        end
    end
    # beam-splitter 
    if type == :bs
        if length(params) == 0
            return sympop(1, 2, (args::AbstractArray)->beamSplitter(args[1]))
        elseif length(params) == 2
            mode1 = params[1]
            mode2 = params[2]
            return sympop(1, max(mode1, mode2),
                          (args::AbstractArray)->beamSplitter(args[1], mode1, mode2))
        elseif length(params) == 3
            mode1 = params[1]
            mode2 = params[2]
            nmodes = params[3]
            return sympop(1, nmodes,
                          (args::AbstractArray)->beamSplitter(args[1], mode1, mode2, nmodes))
        end
    end
    # amplifier
    if type == :amp
        if length(params) == 0
            return sympop(1, 2, (args::AbstractArray)->amplifier(args[1]))
        elseif length(params) == 2
            mode1 = params[1]
            mode2 = params[2]
            return sympop(1, max(mode1, mode2),
                          (args::AbstractArray)->amplifier(args[1], mode1, mode2))
        elseif length(params) == 3
            mode1 = params[1]
            mode2 = params[2]
            nmodes = params[3]
            return sympop(1, nmodes,
                          (args::AbstractArray)->amplifier(args[1], mode1, mode2, nmodes))
        end
    end
    return throw("unimplemented")
end

function (s::SympOp)(args::AbstractArray)
    if length(args) != s.nargs
        return throw(ArgumentError("invalid input"))
    end
    return s.fn(args)
end

function Base.:*(x::SympOp, y::SympOp)
    if x.nmodes != y.nmodes
        return throw(ArgumentError("unequal number of modes"))
    end
    fn = (args::AbstractArray) -> begin
        x.fn(args[1:x.nargs]) * y.fn(args[x.nargs+1:end])
    end
    return sympop(x.nargs + y.nargs, x.nmodes, fn; ensure=x.ensure || y.ensure)
end
function Base.:*(x::SympOp, y::SympOp, z::Vararg{SympOp})
    return reduce(*, [x, y, z...])
end

function Base.transpose(s::SympOp)
    return sympop(s.nargs, s.nmodes, (args::AbstractArray)->transpose(s.fn(args)); ensure=s.ensure)
end

function Base.inv(s::SympOp)
    sympform = sympop(Ω(2*s.nmodes))
    return sympform*transpose(s)*transpose(sympform)
end

function Base.one(s::SympOp)
    return sympop(I(2*s.nmodes))
end
  
# export toSympType, SymplecticOperation
# export SympType, BeamSplitter, PhaseShifting, Custom, Undef

# const Params = AbstractArray{Any}

# @enum SympType PhaseShifting BeamSplitter Custom Undef

# function toSympType(name::AbstractString)
#     if name == "beam_splitter"
#         return BeamSplitter
#     elseif name == "phase_shifting"
#         return PhaseShifting
#     elseif name == "custom"
#         return Custom
#     else
#         return Undef
#     end
# end

# function funcBuild(type::SympType, params::Params)::Func
#     l = length([0 for p in params if p == "_"])

#     func = if type == PhaseShifting
#         phaseShifting
#     elseif type == BeamSplitter
#         beamSplitter
#     elseif type == Custom
#         customReshape
#     else
#         throw(ArgumentError("Operation undefined"))
#     end

#     f = (xs::Vararg) -> begin
#         i = 1
#         ys = []
#         for p in params
#             if p == "_"
#                 push!(ys, xs[i])
#                 i += 1
#             else
#                 push!(ys, p)
#             end
#         end
#         func(ys...)
#     end
#     return Func(l, f)
# end


# struct SymplecticOperation
#     n::Int64
#     Op::Func
#     function SymplecticOperation(n::Int64, Op::Func)
#         return new(n, Op)
#     end
# end

# function (s::SymplecticOperation)(xs::Vararg)
#     return s.Op(xs...)
# end

# function SymplecticOperation(n::Int64, type::SympType, params::Params)::SymplecticOperation
#     op = funcBuild(type, params)
#     return SymplecticOperation(n, op)
# end

# function SymplecticOperation(m::AbstractMatrix)
#     n = size(m, 1) ÷ 2
#     l = 0
#     f = (xs::Vararg) -> m
#     return SymplecticOperation(n, Func(l, f))
# end

# function dsum(s1::SymplecticOperation, s2::SymplecticOperation)::SymplecticOperation
#     n = s1.n + s2.n
#     l1, l2 = s1.Op.l, s2.Op.l
#     l = l1 + l2
#     op = (xs::Vararg) -> begin
#         dsum(s1(xs[1:l1]...), s2(xs[l1+1:l]...))
#     end
#     return SymplecticOperation(n, Func(l, op))
# end

# function dsum(s1::SymplecticOperation, n2::Int64)::SymplecticOperation
#     n = s1.n + n2
#     l = s1.Op.l
#     op = (xs::Vararg) -> begin
#         dsum(s1(xs...), I(2 * n2))
#     end
#     return SymplecticOperation(n, Func(l, op))
# end

# function dsum(n1::Int64, s2::SymplecticOperation)::SymplecticOperation
#     n = n1 + s2.n
#     l = s2.Op.l
#     op = (xs::Vararg) -> begin
#         dsum(I(2 * n1), s2(xs...))
#     end
#     return SymplecticOperation(n, Func(l, op))
# end

# function Base.:*(s1::SymplecticOperation, s2::SymplecticOperation)::SymplecticOperation
#     if s1.n < s2.n
#         s1 = dsum(s1, s2.n - s1.n)
#     end
#     if s1.n > s2.n
#         s2 = dsum(s2, s1.n - s2.n)
#     end

#     n = s1.n
#     op = s1.Op * s2.Op

#     return SymplecticOperation(n, op)
# end

# function Base.inv(s::SymplecticOperation)
#     n = s.n
#     l = s.Op.l
#     op = (xs::Vararg) -> begin
#         -Ω(2 * n) * transpose(s(xs...)) * Ω(2 * n)
#     end
#     return SymplecticOperation(n, Func(l, op))
# end

# function teleportation(s::SymplecticOperation, inModes, outModes)::SymplecticOperation
#     l = s.Op.l
#     n = length(inModes)
#     if n == 0
#         n = s.n
#     end
#     op = (xs::Vararg) -> begin
#         teleportation(s(xs...), inModes, outModes)
#     end
#     return SymplecticOperation(n, Func(l, op))
# end

# function nonSymplecticity(s::SymplecticOperation)::Func
#     l = s.Op.l
#     f = (xs::Vararg) -> begin
#         nonSymplecticity(s(xs...))
#     end
#     return Func(l, f)
# end
