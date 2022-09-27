using LinearAlgebra: I

export circulator

"""
    circulator(perm...)

Return a symplectic matrix permutating modes.

Input can be a vector.
"""
circulator(perm...) = circulator(perm)
circulator(perm::Vector) =
       hvcat(length(perm), [i == perm[j] ? I(2) : zeros(2, 2) for i in 1:length(perm), j in 1:length(perm)]...)
