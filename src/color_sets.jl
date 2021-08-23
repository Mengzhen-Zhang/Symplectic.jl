export sympToGraph,  getColors,  getColorSetsFromSymp, contract

# Generate a Graph based on the given Symplectic Matrix
const Graph = AbstractMatrix{Bool}
const Colors = Vector
const Vertex = Int
const Colored = Bool
const NofModes = Int

function sympToGraph( S::Symp )::Graph
    n = nModes(S)
    return [ norm(S.S[2*i-1:2*i, 2*j-1:2*j]) > tolerance for i in 1:n, j in 1:n ]
end

function colorSuccessors(sucs, Cs::Colors)::Colors
    firstColor = Cs[sucs[1]]
    Cs[sucs] .= firstColor
    return Cs
end  

function  getColors( G::Graph )::Colors
    colored = false
    n = isSquare(G) # = length(Cs)
    Cs = [1:n;]
    for v in 1:n
        sameColorVs = G[ findall(x->x == Cs[v], Cs ), :]
        sucs = findall(vi->any(sameColorVs[:, vi] ),  [1:n;] )
        nColors = length( Set(Cs[sucs]) )
        if nColors > 1
            Cs = colorSuccessors( sucs, Cs )
            colored = true
        end
    end
    return colored ? getColors(G) : Cs
end

# Is i a successor of j (both colored)?
function areColoredVeticesConnected(ci, cj, colorSets::Dict, G::Graph)::Bool
    vis = colorSets[ci]
    vjs = colorSets[cj]
    return any(G[vis, vjs])
end

getAllColors(Cs::Colors)::Vector = [Set(Cs)...]

function getColorSets(Cs::Colors)::Dict
    allColors = getAllColors(CS)
    return Dict(c => findall(x -> x == c, Cs) for c in allColors)
end

getColorSetsFromSymp = getColorSets ∘ getColors ∘ sympToGraph

function contract( G::Graph, Cs::Colors )
    colorSets = getColorSets(Cs)
    return [areColoredVeticesConnected(ci, cj, colorSets, G) for  ci in allColors, cj in allColors]
end
