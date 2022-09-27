module Symplectic

include("SymplecticForm.jl")
include("inverse_symplectic_cayley_transform.jl") # <- SymplecticForm
include("symplectic_cayley_transform.jl")         # <- SymplecticForm
include("non_symplecticity.jl")                   # <- SymplecticForm

include("dsum.jl")                             
include("phase_shifting.jl")                    # <- dsum
include("amplifier.jl")                         # <- dsum
include("extend.jl")                            # <- dsum
include("beam_splitter.jl")                     # <- dsum, phase_shifting
include("random_phase_shifting.jl")             # <- phase_shifting
include("local_symplectic_matrix.jl")           # <- dsum, phase_shifting
include("random_local_symplectic_matrix.jl")    # <- local_symplectic_matrix

include("circulator.jl")                        

include("feedforward.jl")

include("teleportation.jl")

include("interference_based_sequence.jl")

include("to_qpqp_basis.jl")             
include("to_qqpp_basis.jl")
include("pre_iwasawa_factorization.jl") # <- to_qqpp_basis, to_qpqp_basis
include("spectrum.jl")                  # <- to_qqpp_basis, to_qpqp_basis

include("helper.jl")
include("Func.jl")
include("SymplecticOperation.jl")
include("SymplecticCircuit.jl")

export SymplecticForm, Ω,
       SymplecticOperation,
       toSympType
end
