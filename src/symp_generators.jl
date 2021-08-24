export Id, monoSymp, beamSplitter, circulator
export randomSymmetric, randomSymp, randomGenericSymp
export randomPassiveLocalSymp, randomLocalSymp, randomizeSymp

# From Number of Modes to Identity Symplectic Matrix
Id(n::Int)::Symp = Matrix(1I, 2*n, 2*n)

monoSymp(θ::Real)::Symp =  exp(θ*Ω(1)) 
# Generte a Two-Mode BeamSplitter
beamSplitter(angle::Real)::Symp = monoSymp(angle) ⊗ Id(1)
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
randomPassiveLocalSymp(n::NofModes)::Symp = n == 1 ? monoSymp(2π*rand()) : ⊕([monoSymp(2π*rand()) for i in 1:n]...)

randomLocalSymp(n::NofModes, range::Real = 3)::Symp = n == 1 ? randomSymp(1, range) : ⊕([randomSymp(1, range) for i in 1:n]...)

function randomizeSymp(S::Symp, l::Int)::Symp
    n = nModes(S)
    if l == 1
        return randomPassiveLocalSymp(n) * S  * randomPassiveLocalSymp(n)
    else
        return *([randomizeSymp(S, 1) for i in 1:l]...)
    end
end


