# depend on phase_shfiting.jl
# depend on dsum.jl

export localSymplecticMatrix

localSymplecticMatrix(θ::Real, l::Real, ϕ::Real) =
    phaseShifting(θ) * [l 0; 0 1/l] * phaseShifting(ϕ)

"""
    localSymplecticMatrix(θs::Union{Real,Vector}, ls::Union{Real,Vector}, ϕs::Union{Real,Vector})

Return a direct sum of 2×2 `localSymplecticMatrix`s according to Euler decompostion.

Return a 2×2 matrix when the input arguments are real numbers.
"""
localSymplecticMatrix(θs::Vector, ls::Vector, ϕs::Vector) =
    dsum(map(args -> localSymplecticMatrix(args...), zip(θs, ls, ϕs))...)
