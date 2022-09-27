# depend on phase_shfiting.jl

export randomPhaseShifiting

"""
    randomPhaseShifiting(n::Int)

Return an n-mode random phase-shifting symplectic matrix.
"""
randomPhaseShifiting(n::Int) =
       n == 1 ? phaseShifting(2π * rand()) : phaseShifting((2π * rand(n))...)
