module Symplectic
    import LinearAlgebra
    
    include("Utilities/LinearAlgebraUtilities.jl")
    import .LinearAlgebraUtilities

    include("DataTypes/SymplecticMatrix.jl")

    include("DataTypes/MatrixSequence.jl")



    # Data type: SymplecticMatrix
    export SymplecticMatrix,
           Omega, Ω,
           isGeneric,
           ⊗,
           dsum, ⊕

    # Data type: MatrixSequence
    export MatrixSequence,
           sliceInner,
           slice,
           setInterspersion,
           setSequence,
           replaceSequence

    # Utility module for linear algebra
    export LinearAlgebraUtilities

end