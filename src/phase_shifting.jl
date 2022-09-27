# depend on dsum.jl

export phaseShifting

"""
    phaseShifting(θ...)

Return a (multi-)mode phase-shifting symplectic matrix.

Size of the returned matrix is specified by the length of the input.
"""
phaseShifting(θ::Real) = [cos(θ) sin(θ); -sin(θ) cos(θ)]
phaseShifting(θs...) = dsum(map(phaseShifting, θs)...)
