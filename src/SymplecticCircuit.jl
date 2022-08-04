import JSON
import Optim: optimize, minimizer
import ReverseDiff: gradient

function buildCircuit(dict)
    n = dict["number_of_modes"]
    ops = dict["symplectic_operations"]

    if length(ops) < 1
        throw("Operation list is empty")
    end

    nameToFunction = Dict(
        "beam_splitter" => (θ, i, j) -> beamSplitter(θ, i, j, n),
        "phase_shifting" =>  phaseShifting,
        "custom" => (params...) -> reshape([params...], (2*n, 2*n))
    )

    sympList = []
    l = 0

    for op in ops
        name = op["name"]
        if name == "_"
            l += (2*n)^2
        else
            params = op["params"]
            if name == "beam_splitter" && params[1] == "_"
                l += 1
            elseif name == "phase_shifting" && params[1] == "_"
                l += n
            end
        end
    end


    f = (argVec) -> begin
            i = 1
            output = I(2*n)
    
            for op in ops
                name = op["name"]
                if name == "_"
                    newOp = reshape(argVec[i:i+(2*n)^2-1], (2*n, 2*n))
                    sympList = vcat(sympList, [i])
                    i += (2*n)^2
                else
                    params = op["params"]
                    if name == "beam_splitter" && params[1] == "_"
                        newOp = beamSplitter(argVec[i], params[2], params[3], n)
                        i += 1
                    elseif name == "phase_shifting" && params[1] == "_"
                        newOp = phaseShifting(argVec[i:i+n-1]...)
                        i += n
                    else
                        newOp = nameToFunction[name](params...)
                    end
                end
                output = newOp * output
            end    
            l = length(argVec)
            return output
        end

    return f, n, l, sympList
end

function buildCircuitFromJSON(json)
    dict = JSON.parse(json)
    return buildCircuit(dict)
end

function buildCircuitFromFile(filename)
    dict = JSON.parsefile(filename)
    return buildCircuit(dict)
end

function buildCostFromFile(filename)
    f, n, l, lst = buildCircuitFromFile(filename)
    cost = (xs) -> begin
        output = (x -> tr(transpose(x) * x))(f(xs))
        for s in lst
            output += (x -> tr(transpose(x) * x))(reshape(xs[s:s+(2*n)^2-1], (2*n, 2*n)))
        end
        return output
    end
    return cost, l
end


f, n, l, lst = buildCircuitFromFile("test/test.json")

function findSolutionFromFile(filename)
    cost, l = buildCostFromFile("test/test.json")
    function df!(G, xs)
        Gs = gradient(cost, xs)
        for i in 1:length(xs)
            G[i] = Gs[i]
        end
    end

    res = optimize(cost, df!,
                    rand(l),
                    method = BFGS(),
                    iterations = 1000)

    return res
end

    
res = findSolutionFromFile("test/test.json")
minimizer(res)