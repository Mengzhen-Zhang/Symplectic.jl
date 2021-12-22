# Importing functions from, Adding functions to, and renaming functions in LinearAlgebra

module LinearAlgebraUtilities
    import LinearAlgebra.norm,
           LinearAlgebra.checksquare

    export norm,
           isSquare,
           isGeneric,
           isSymplectic,
           Omega, Ω,
           dsum, ⊕

    Omega(n::Int)::AbstractMatrix = cat(fill([0 1; -1 0], n)...; dims=(1,2)) 
    Omega(matrix::AbstractMatrix)::AbstractMatrix = Omega(size(matrix)[1] ÷ 2)
    Ω = Omega
    
    function isSquare(matrix::AbstractMatrix)::Bool
        try 
            n = checksquare(matrix)
            return true
        catch err
            return false
        end
    end
    isGeneric(matrix::AbstractMatrix, tolerance::Real = 1e-6)::Bool = all(x -> abs(x) > tolerance, matrix )
    function isSymplectic(matrix::AbstractMatrix, tolerance = 1e-6::Real)::Bool
        if !(isSquare(matrix)) return false end
        if (checksquare(matrix) % 2 != 0) return false end
        if (norm(matrix*Omega(checksquare(matrix) ÷ 2)*matrix' - Omega(checksquare(matrix) ÷ 2)) >= tolerance) return false end
        return true
    end

    dsum(M1::AbstractMatrix, M2::AbstractMatrix)::AbstractMatrix = cat(M1, M2; dims=(1,2))
    dsum(M1::AbstractMatrix, M2::AbstractMatrix, Ms::Vararg{AbstractMatrix})::AbstractMatrix = dsum(dsum(M1, M2),  mapreduce(identity, dsum, Ms ))
    ⊕ = dsum
end