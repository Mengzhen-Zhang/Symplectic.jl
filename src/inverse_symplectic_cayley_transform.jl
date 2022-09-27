# depend on SymplecticForm.jl

using LinearAlgebra: I

export inverseSymplecticCayleyTransform, isct

"""
    inverseSymplecticCayleyTransform(S::AbstractMatrix)

return `-Ω*((I-S)^-1 - I/2)`.
"""
inverseSymplecticCayleyTransform(S::AbstractMatrix) = -Ω * ((I - S)^-1 - I / 2)
const isct = inverseSymplecticCayleyTransform
