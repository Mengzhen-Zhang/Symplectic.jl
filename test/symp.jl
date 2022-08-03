using Symplectic
using LinearAlgebra
using ForwardDiff
using Optim

const bs = beamSplitter
const ϕ = π/4
θ = π / 4
invBS = inv(beamSplitter(π/4)) ⊕ I(2*2)
S = bs(ϕ, 1, 2, 4)*bs(θ, 2, 3, 4)*bs(ϕ, 3, 4, 4)*bs(ϕ, 1, 4, 4)*bs(ϕ, 1, 2, 4)*bs(θ, 2, 3, 4)*bs(ϕ, 3, 4, 4)*bs(ϕ, 1, 4, 4)
S = *(generify(S; randomize = randomPhaseShifiting )...)

dseq2 = decoupleSequence([S for i in 1:4]...);
dseq1 = decoupleSequence([i==4 ? S*invBS : S for i in 1:4]...);
S2 = *(dseq2...);
S1 = *(dseq1...);
S2 = S2[1:2, 1:2] ⊕ *(generify( S2[3:end, 3:end], randomize = randomLocalSymplecticMatrix)...)



genericity(S2[3:end, 3:end])
round.(S2[3:end, 3:end]; digits = 3)

seq = ibs(S; T = BS);
round.(*(seq...); digits = 3)
round.(seq[5]; digits=3)

tmp = [ 1 1 1 1;
        -1 1 -1 1;
        -1  -1 1 1;
        1 -1 -1 1] / 2

nonSymplecticity(tmp)                

transpose(tmp) * Ω * tmp

a, b, c = preIwasawaFactorization(tmp)
a
b
c

lstf(tmp[:,1], tmp[1,:])

S * lstf(S[:, 1], S[1, :]) * S

opnorm(round.(transpose(S) * S - I))


function func(θs::Vector)
        T = teleportation(phaseShifting([0, 0, θs[3:4]...]...)*Matrix(S)*phaseShifting([θs[1:2]..., 0, 0]...), [1, 2], [1, 2])
        return tr((T - I(4)) * transpose(T - I(4)))
end

function func2(θs::Vector)
        teleportation(phaseShifting([0, 0, θs[3:4]...]...)*Matrix(S)*phaseShifting([θs[1:2]..., 0, 0]...), [1, 2], [1, 2])
end



g = x -> ForwardDiff.gradient(func, x)

function g2(x::Vector, storage::Vector)
        tmp = g(x)
        for i in 1:length(storage)
                storage[i] = tmp[i]
        end
end

g([1,2,0,0])

func([1,2,0,0])

res = optimize(func, g2, [1.0, 2.0, 3.0, 4.0], SimulatedAnnealing())
func2(Optim.minimizer(res))