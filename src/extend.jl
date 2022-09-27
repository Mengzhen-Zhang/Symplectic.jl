# depend on dsum.jl

using LinearAlgebra: I

export extend

"""
    extend(m::AbstractMatrix, pos::AbstractVector{Int}, modes::Int)

Return a direct sum, `A`, of `m` and an identity matrix, with the `pos[i]`th mode of `A` equal to
the `i`th mode of the input matrix `m`.
"""
function extend(m::AbstractMatrix, pos::AbstractVector{Int}, modes::Int)
    lst = zeros(Int64, modes)
    for (i, p) in enumerate(pos)
        lst[p] = i
    end
    j = size(m, 1) + 1
    for i in 1:modes
        if lst[i] == 0
            lst[i] = j
            j += 1
        end
    end
    return dsum(m, I(modes - size(m, 1)))[lst, lst]
end
