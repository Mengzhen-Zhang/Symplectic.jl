export sympCayleyTransform, inverseSympCayleyTransform, sympRound

# Symplectic Cayley Transform: Output a symplectic matrix when the input is a symmetric matrix
function sympCayleyTransform(M::AbstractMatrix)::AbstractMatrix
    n = isSquare(M) ÷ 2
    Ω = cat(fill([0 1; -1 0], n)...; dims=(1,2))
    return I - (Ω * M  + I / 2  )^-1
end
# Invese Symplectic Cayley Transform: Output a symmetric matrix when the input is a symplectic matrix
function inverseSympCayleyTransform(S::AbstractMatrix)::AbstractMatrix
    n = isSquare(S) ÷ 2
    Ω = cat(fill([0 1; -1 0], n)...; dims=(1,2))
    return Symmetric(-Ω*((I - S)^-1 - I/2))
end
# Rounding a Matrix using Cayley Transform
function sympRound(S::AbstractMatrix)::AbstractMatrix
    n = isSquare(S) ÷ 2
    Ω = cat(fill([0 1; -1 0], n)...; dims=(1,2))
    if det(S - I) < tolerance
        return -Ω * ( (sympCayleyTransform ∘ inverseSympCayleyTransform)(Ω * S) )
    else
        return (sympCayleyTransform ∘ inverseSympCayleyTransform)( S )
    end
end

# Rounding a Symplectic Matrix using Cayley Transform
sympRound(S::Symp)::Symp = sympRound(S.S)