# depend on local_symplectic_matrix.jl

export randomLocalSymplecticMatrix

"""
    randomLocalSymplecticMatrix(n::Int, bound::Real=0.2)

Return a random direct sum of 2×2 `localSymplecticMatrix`s.

`bound` bounds the degree of squeezing.
"""
randomLocalSymplecticMatrix(n::Int; bound::Real=0.2) =
       n == 1 ? localSymplecticMatrix(2π * rand(), bound * (rand() - 1 / 2) + 1.0, 2π * rand()) :
       localSymplecticMatrix(2π * rand(n), bound * (rand(n) .- 1 / 2) .+ 1.0, 2π * rand(n))
