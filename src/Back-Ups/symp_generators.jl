export Id, monoSymp, beamSplitter, circulator
export randomSymmetric, randomSymp, randomGenericSymp
export randomPassiveLocalSymp, randomLocalSymp, randomizeSymp

# From Number of Modes to Identity Symplectic Matrix
Id(n::Int)::Symp = Matrix(1I, 2*n, 2*n)

monoPassiveSymp(θ::Real)::Symp =  exp(θ*Ω(1))
monoSymp(θ1::Real, l::Real, θ2::Real)::Symp =  monoPassiveSymp(θ1) * Symp([l 0 ; 0  1/l ]) * monoPassiveSymp(θ2)

# Generte a Two-Mode BeamSplitter
beamSplitter(angle::Real)::Symp = monoPassiveSymp(angle) ⊗ Id(1)
# Generate a BeamSplitter between two modes in a multi-mode System
function beamSplitter(angle::Real, mode1, mode2, nModes)::Symp
    mode1, mode2 = sort([mode1, mode2])
    lst = [ [[2 * (i + 2) - 1, 2 * (i + 2)] for i in 1:mode1-1]... ; [1, 2] ; [[2 * (i + mode1 + 1) - 1, 2 * (i + mode1 + 1)] for i in 1: mode2 - mode1 - 1]...; [3, 4] ; [[2 * (i + mode2 ) - 1, 2 * (i + mode2 )]  for i in 1:nModes-mode2]...  ];
    return dsum( beamSplitter(angle), Id(nModes - 2) )[lst, lst]
end

# Generate a Circulator-Like Symplectic Matrix according to the given permutation
circulator(perm::Vector)::Symp = hvcat(length(perm), [i == perm[j] ? Id(1).S : zeros(2, 2)   for i in 1:length(perm), j in 1:length(perm) ]...)


# From Number of Modes to Random Symmetric Matrix
function randomSymmetric(n::Int, range::Real = 3)::Matrix{BigFloat}
    if range <= 0
        return error("range must be positive")
    end
    return Symmetric( range*(rand(2*n, 2*n) - I ) )
end

# From Number of Modes to Random Symplectic Matrix
function randomSymp(n::Int, range::Real = 3)::Symp
    S1 = Id(n) - (Ω(n) * randomSymmetric(n, range) + Id(n) / 2  )^-1
    S2 = Id(n) - (Ω(n) * randomSymmetric(n, range) + Id(n) / 2  )^-1
    return  S1 * S2 
end
function randomGenericSymp(n::Int, range::Real = 3)::Symp
        m = randomSymp(n, range)
        return isGeneric(m) ? m : randomGenericSymp(n, range)
end
randomPassiveLocalSymp(n::NofModes)::Symp = n == 1 ? monoPassiveSymp(2π*rand()) : ⊕([monoPassiveSymp(2π*rand()) for i in 1:n]...)

randomLocalSymp(n::NofModes, range::Real = 0.2)::Symp = n == 1 ? monoSymp(2π*rand(), range*( rand() -1/2) + 1.0 , 2π*rand()) : ⊕([ monoSymp(2π*rand(), range*( rand() -1/2) + 1.0 , 2π*rand()) for i in 1:n]...)

function randomizeSympPassive(S::Symp, l::Int)::Symp
    n = nModes(S)
    if l == 1
        return randomPassiveLocalSymp(n) * S  * randomPassiveLocalSymp(n)
        # return randomLocalSymp(n) * S  * randomLocalSymp(n)
    else
        return *([randomizeSympPassive(S, 1) for i in 1:l]...)
    end
end

function randomizeSymp(S::Symp, l::Int)::Symp
    n = nModes(S)
    if l == 1
        return randomLocalSymp(n) * S  * randomLocalSymp(n)
        # return randomLocalSymp(n) * S  * randomLocalSymp(n)
    else
        return *([randomizeSymp(S, 1) for i in 1:l]...)
    end
end

