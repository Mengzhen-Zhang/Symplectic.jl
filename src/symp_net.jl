# depend on symp_op.jl

struct SympStack
    symp_ops::AbstractArray{SympOp}
    function SympStack(symp_ops::AbstractArray{SympOp})
        if length(symp_ops) == 0
            throw(ArgumentError("symp_ops is empty"))
        end
        for i in 1:length(symp_ops)
            if symp_ops[i].nmodes != symp_ops[1].nmodes
                throw(ArgumentError("unequal number of modes"))
            end
        end
        return new(symp_ops)
    end
end

sympstack(symp_ops::AbstractArray{SympOp}) = SympStack(symp_ops)

Base.convert(::Type{AbstractArray{SympOp}}, net::SympStack) = net.symp_ops

Base.length(net::SympStack) = length(net.symp_ops)

Base.size(net::SympStack) = (length(net), )

function Base.getindex(net::SympStack, i::Int)
    if i < 1 || i > length(net)
        throw(BoundsError(net, i))
    end
    return net.symp_ops[i]
end
function Base.getindex(net::SympStack, n::AbstractRange{<:Integer})
    v = []
    @inbounds for (i, ii) in enumerate(n)
        push!(v, ii)
    end
    return v
end

function Base.setindex!(net::SympStack, s::SympOp, i::Int)
    if i < 1 || i > length(net)
        throw(BoundsError(net, i))
    end
    if s.nmodes != net[1].nmodes
        throw(ArgumentError("unequal number of modes"))
    end
    @inbounds net.symp_ops[i] = s
    return net[i]
end
function Base.setindex!(net::SympStack, X::Union{SympStack, AbstractArray{SympOp}}, n::AbstractRange{<:Integer})
    @inbounds for (i, ii) in enumerate(n)
        if net[1].nmodes != X[i].nmodes
            throw(ArgumentError("unequal number of modes"))
        end
        net[ii] = X[i]
    end
    return net[n]
end

Base.firstindex(net::SympStack) = 1

Base.lastindex(net::SympStack) = length(net)

function Base.push!(net::SympStack, s::SympOp)
    push!(net.symp_ops, s)
    return net
end
function Base.push!(net::SympStack, s::Vararg{SympOp})
    push!(net.symp_ops, s...)
    return net
end

function Base.append!(net::SympStack, X::SympStack)
    append!(net.symp_ops, X.symp_ops)
    return net
end
function Base.append!(net::SympStack, X::AbstractArray{SympOp})
    append!(net.symp_ops, X)
    return net
end
function Base.append!(net::SympStack, X::Vararg{Union{SympStack, AbstractArray{SympOp}}})
    for x in X
        append!(net, x)
    end
    return net
end

Base.vcat(net::SympStack) = net
function Base.vcat(net::Union{SympStack, AbstractArray{SympOp}}, X::Vararg{Union{SympStack, AbstractArray{SympOp}}})
    flag = false
    if isa(net, SympStack)
        flag = true
    end
    v = [convert(AbstractArray{SympOp}, net)]
    for x in X
        if !flag && isa(x, SympStack)
            flag = true
        end
        push!(convert(AbstractArray{SympOp}, x))
    end
    return sympstack(v)
end
