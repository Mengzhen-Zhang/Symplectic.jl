import JSON
import Optim
import ReverseDiff: gradient
import LinearAlgebra: tr

export SymplecticCircuit, buildCircuit, buildCircuitFromFile, buildCircuitFromJSON, buildConstraint, buildConstraintFromFile, findSolutionFromFile

struct SymplecticCircuit
    circuit::AbstractArray{SymplecticOperation}
    inModes::AbstractArray
    outModes::AbstractArray
end

function regularize(n::Int64, type::SympType, params::Params)::Params
    if type == PhaseShifting
        return params
    elseif type == BeamSplitter
        return Any[params..., n]
    elseif type == Custom
        return Any[params..., 2*n]
    else
        throw(ArgumentError("Operation undefined"))
    end
end


function getArguments(n::Int64, op::Dict)
    type = toSympType(op["name"])
    params = get(op, "params",
        if type == Custom
            params = Any["_" for i in 1:2*n*2*n]
        elseif type == PhaseShifting
            params = Any["_" for i in 1:n]
        end)
    return n, type, regularize(n, type, params)
end

function buildCircuit(dict::Dict)::SymplecticCircuit
    n = get(dict, "number_of_modes", -1)
    if n == -1
        throw("must specify number_of_modes")
    end

    ops = get(dict, "symplectic_operations", [])
    if length(ops) == 0
        throw("operation list cannot be empty")
    end

    inModes = get(dict, "input_modes", [])
    outModes = get(dict, "output_modes", [])

    circuit = [SymplecticOperation(getArguments(n, op)...) for op in ops[end:-1:1]]

    return SymplecticCircuit(circuit, inModes, outModes)
end

function buildCircuitFromJSON(json::AbstractString)::SymplecticCircuit
    dict = JSON.parse(json)
    return buildCircuit(dict)
end

function buildCircuitFromFile(filename::AbstractString)::SymplecticCircuit
    dict = JSON.parsefile(filename)
    return buildCircuit(dict)
end

function buildConstraint(sc::SymplecticCircuit, target=nonSymplecticity)::CustomFunction
    so = reduce(*, sc.circuit)
    l = so.Op.l
    f = (xs::Vararg) -> begin
        constraint = nonSymplecticity(so(xs...))
        lcur = 1
        for op in sc.circuit
            lop = op.Op.l
            if lop > 0
                constraint += nonSymplecticity(op(xs[lcur:lcur+lop-1]...))
            end
            lcur += lop
        end
        tel = teleportation(so(xs...), sc.inModes, sc.outModes)
        constraint += target(tel)
        
        return constraint
    end

    return CustomFunction(l, f)
end

function buildConstraintFromFile(filename::AbstractString, target=nonSymplecticity)::CustomFunction
    sc = buildCircuitFromFile(filename)
    return buildConstraint(sc, target)
end

function findSolutionFromFile(filename, target=nonSymplecticity)
    constraint = buildConstraintFromFile(filename, target)

    l = constraint.l

    f(xs) = constraint(xs...)

    function df!(G, xs)
        grad = gradient(f, xs)
        for i in 1:length(xs)
            G[i] = grad[i]
        end
    end

    res = Optim.optimize(f, df!,
        rand(l),
        method=Optim.LBFGS(),
        iterations=1000)

    return res
end