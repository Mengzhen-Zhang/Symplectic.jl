abstract type AbstractFunc end

struct Func
    l::Integer
    f::Function
end

convert(Func, fn::AbstractFunc) = Func(get_l(fn), get_f(fn))

get_l(fn::AbstractFunc) = get_l(fn::Func) # Interface
get_f(fn::AbstractFunc) = get_f(fn::Func) # Interface

get_l(fn::Func) = fn.l
get_f(fn::Func) = fn.f


Base.:+(fn::AbstractFunc) = fn

Base.:-(fn1::AbstractFunc, fn2::AbstractFunc) = fn1 + (- fn2)

function Base.:-(fn::Func)
    l, f = get_l(fn), get_f(fn)
    minus_f = (xs::Vararg) -> - f(xs...)
    return Func(l, minus_f)
end

function Base.:+(fn1::Func, fn2::Func)
    l1, f1 = get_l(fn1), get_f(fn1)
    l2, f2 = get_l(fn2), get_f(fn2)
    l = l1 + l2
    sum_f = (xs::Vararg) -> begin
	x1s, x2s = xs[1:l1], xs[l1+1:l]
        return f1(x1s...) + f2(x2s...)
    end
    return Func(l, sum_f)
end

function Base.:*(fn1::Func, fn2::Func)
    l1, f1 = get_l(fn1), get_f(fn1)
    l2, f2 = get_l(fn2), get_f(fn2)
    l = l1 + l2
    sum_f = (xs::Vararg) -> begin
	x1s, x2s = xs[1:l1], xs[l1+1:l]
        return f1(x1s...) * f2(x2s...)
    end
    return Func(l, sum_f)
end

function dsum(fn1::Func, fn2::Func)
    l1, f1 = get_l(fn1), get_f(fn1)
    l2, f2 = get_l(fn2), get_f(fn2)
    l = l1 + l2
    dsummed_f = (xs::Vararg) -> begin
	x1s, x2s = xs[1:l1], xs[l1+1:l]
        return dsum(f1(x1s...), f2(x2s...))
    end
    return Func(l, dsummed_f)
end

function dsum(fn1, fn2::Func)
    l1 = 0
    l2, f2 = get_l(fn2), get_f(fn2)
    l = l1 + l2
    dsummed_f = (xs::Vararg) -> begin
	x1s, x2s = xs[1:l1], xs[l1+1:l]
        return dsum(fn1, f2(x2s...))
    end
    return Func(l, dsummed_f)
end

function dsum(fn1::Func, fn2)
    l1, f1 = get_l(fn1), get_f(fn1)  
    l2 = 0
    l = l1 + l2
    dsummed_f = (xs::Vararg) -> begin
	x1s, x2s = xs[1:l1], xs[l1+1:l]
        return dsum(f1(x1s...), fn2)
    end
    return Func(l, dsummed_f)
end



function Base.transpose(fn::Func)
    l, f = get_l(fn), get_f(fn)
    transposed_f = (xs::Vararg) -> transpose(f(xs...))
    return Func(l, transposed_f)
end

function (fn::AbstractFunc)(xs::Vararg)
    l, f = get_l(fn), get_f(fn)
    if length(xs) < l
        throw(ArgumentError("not enough arguments"))
    end
    return f(xs...)
end
