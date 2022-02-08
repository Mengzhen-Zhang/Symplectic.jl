using Symplectic
using LinearAlgebra
using ForwardDiff
using Optim

const bs = beamSplitter
const ϕ = π/4
const θ = π / 4
# S = bs(ϕ, 1, 2, 5)*bs(ϕ, 2, 3, 5)*bs(θ, 3, 4, 5)*bs(ϕ, 4, 5, 5)*bs(ϕ, 1, 5, 5)
# S = bs(ϕ, 1, 2, 7)*bs(ϕ, 2, 3, 7)*bs(ϕ, 3, 4, 7)*bs(ϕ, 4, 5, 7)*bs(θ, 5, 6, 7)*bs(ϕ, 6, 7, 7)*bs(ϕ, 1, 7, 7)
# S = S * S * S
# S = sct(Symmetric(0.2*( rand(14, 14) .-1/2) .+ 1.0))

function evaluation(θs::Vector)
    S = solution(θs)
    T = [0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0] # bs(π/2) # I(4) # bs(π/3) * Diagonal([2, 1/2, 1/3, 3]) * bs(π/3)
    return sqrt(tr((S - T) * transpose(S - T)))
end

function solution(θs::Vector)
    return teleportation(phaseShifting([1, 1, θs[1:5]...]...)*Matrix(S)*phaseShifting([1, 1, θs[6:10]...]...), [1, 2], [1, 2])
end

res = optimize(evaluation,
                2*π*rand(10),
                method = BFGS(), 
                iterations = 1000000,
                autodiff = :forward)
round.(solution(Optim.minimizer(res)), digits = 4)