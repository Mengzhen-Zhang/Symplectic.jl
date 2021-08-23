#= 
    Install Package:  
        In Julia REPL:
            1. type  ']'  to enter the  'pkg'  interface
            2. type  ' add "https://github.com/Mengzhen-Zhang/Symplectic.jl.git" '
=#
module Symplectic

using LinearAlgebra

const tolerance = 10^-6

const Graph = AbstractMatrix{Bool}
const Colors = Vector
const Vertex = Int
const Colored = Bool
const NofModes = Int

# Define DataType Symp for Symplectic Matrices
include("symp.jl")

# Functions for Generating Symplectic Matrices
include("symp_generators.jl")

# Functions for Generating Color Sets
include("color_sets.jl")

# Implemention of Interfernce Based Algorithms
include("interference_based_sequence.jl")

# Implemention of Teleportation Based Symplectic Matrix Generator
include("teleportation_based_generator.jl")

# Implemention of (Inverse) Symplectic Cayley Transform
include("cayley_transform.jl")

end
