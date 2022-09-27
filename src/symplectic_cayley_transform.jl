# depend on SymplecticForm.jl

using LinearAlgebra: I

export symplecticCayleyTransform, sct

"""
    symplecticCayleyTransform(M::AbstractMatrix)

Return `I - (Ω * M + I / 2)^-1`.
"""
symplecticCayleyTransform(M::AbstractMatrix) = I - (Ω * M + I / 2)^-1
const sct = symplecticCayleyTransform
