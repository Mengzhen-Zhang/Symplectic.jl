module Symplectic

using LinearAlgebra:
       checksquare, opnorm, eigen, isposdef, issymmetric, norm, dot

include("SymplecticForm.jl")
include("helper.jl")
include("SymplecticOperation.jl")
include("SymplecticCircuit.jl")
include("CoupledResonators.jl")

export SymplecticForm, Ω,

       # interfaces
       checkSymplectic, modes, getQ, getP, ⊗, dsum, ⊕,

       # constructors
       localSymplecticTransformation, lstf,
       randomPhaseShifiting,
       localSymplecticMatrix,
       randomLocalSymplecticMatrix,
       beamSplitter,
       circulator,
       generify,
       phaseShifting,

       # algorithms and formulae implementation
       symplecticCayleyTransform, sct,
       inverseSymplecticCayleyTransform, isct,
       interferenceBasedSequence, ibs,
       decoupleSequence,
       teleportation,
       preIwasawaFactorization, pif,
       spectrum,

       # helpers
       genericity,
       nonSymplecticity,
       toQPQPBasis,
       toQQPPBasis,

       # symplectic operation
       SymplecticOperation,
       toSympType,

       # generating scatteringMatrix using coupled-mode theory
       CoupledResonators,
       addActive,
       addPassive,
       addGammaEx,
       addGammaIn,
       scatteringMatrix

end