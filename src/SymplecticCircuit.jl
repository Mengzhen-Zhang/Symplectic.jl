import JSON
import Optim: optimize, minimizer, BFGS
import ReverseDiff: gradient
import LinearAlgebra: tr

struct SymplecticCircuit
    circuit ::AbstractArray{SymplecticOperation}
    teleported ::Bool
    inModes::AbstractArray
    outModes::AbstractArray
end 


function regularize(n::Int64, type::SympType, params::AbstractArray)
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

function buildCircuit(dict::Dict) ::SymplecticCircuit
    try
        n = dict["number_of_modes"]    
    catch
        throw("must specify number_of_modes")
    end

    try
        ops = dict["symplectic_operations"]
    catch
        throw("operation list cannot be empty")
    end

    if length(ops) < 1
        throw("operation list cannot be empty")
    end

    try
        inModes = dict["input_modes"]
        outModes = dict["output_modes"]
        teleported = true
    catch
        inModes = []
        outModes = []
        teleported = false
    end

    op = ops[end]
    type = toSympType(op["name"])
    try
        params = regularize(n, type, op["params"])
    catch
        if type == Custom
            params = Any["_" for i in 1:2*n*2*n]
        elseif type == PhaseShifting
            params = Any["_" for i in 1:n]
        end
    end
    circuit = [SymplecticOperation(n, type, Params)]
    
    for op in ops[end-1:-1:1]
        type = toSympType(op["name"])
        try
            params = regularize(n, type, op["params"])
        catch
            if type == Custom
                params = Any["_" for i in 1:2*n*2*n]
            elseif type == PhaseShifting
                params = Any["_" for i in 1:n]
            end
        end 
        push!(circuit, SymplecticOperation(n, type, Params))
    end

    return SymplecticCircuit(circuit, teleported, inModes, outModes)
end

function buildCircuitFromJSON(json::AbstractString) ::SymplecticCircuit
    dict = JSON.parse(json)
    return buildCircuit(dict)
end

function buildCircuitFromFile(filename::AbstractString) ::SymplecticCircuit
    dict = JSON.parsefile(filename)
    return buildCircuit(dict)
end

function buildConstraint(sc::SymplecticCircuit, target=nothing) ::CustomFunction
    so = reduce(*, sc.circuit)
    constraint =  reduce(+, map(nonSymplecticity, sc.circuit))
    constraint += nonSymplecticity(so)
    if sc.teleported
        constraint += nonSymplecticity(teleportation(so, sc.inModes, sc.outModes) - target)
    end
    return constraint
end

function buildConstraintFromFile(filename) ::CustomFunction
    sc = buildCircuitFromFile(filename)
    return buildConstraintFromFile(sc)
end

function findSolutionFromFile(filename)
    constraint = buildConstraintFromFile(filename)

    l, f = constraint.l, constraint.f
    function df!(G, xs)
        Gs = gradient(f, xs)
        for i in 1:length(xs)
            G[i] = Gs[i]
        end
    end

    res = optimize(f, df!,
                    rand(l),
                    method = BFGS(),
                    iterations = 1000)

    return res
end

    
res = findSolutionFromFile("test/test.json")
minimizer(res)