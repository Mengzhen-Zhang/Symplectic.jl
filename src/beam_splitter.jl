# depend on phase_shifting.jl
# depend on dsum.jl

using LinearAlgebra: I

export beamSplitter

"""
    beamSplitter(angle::Real[, mode1, mode2[, modes]])

Return a beam-splitter symplectic matrix, acting on `mode1` and `mode2` of
a system containg `modes` of modes.

Return default to return a 4×4 matrix. `modes` default to `max(mode1, mode2)`.
"""
beamSplitter(angle::Real) = kron(phaseShifting(angle), I(2))
# Generate a BeamSplitter between two modes in a multi-mode System
function beamSplitter(angle::Real, mode1, mode2, modes)
       mode1, mode2 = sort([mode1, mode2])
       lst = [[[2 * (i + 2) - 1, 2 * (i + 2)] for i in 1:mode1-1]...; [1, 2]; [[2 * (i + mode1 + 1) - 1, 2 * (i + mode1 + 1)] for i in 1:mode2-mode1-1]...; [3, 4]; [[2 * (i + mode2) - 1, 2 * (i + mode2)] for i in 1:modes-mode2]...]
       return dsum(beamSplitter(angle), I(2 * modes - 4))[lst, lst]
end
beamSplitter(angle::Real, mode1, mode2) = beamSplitter(angle, mode1, mode2, max(mode1, mode2))
