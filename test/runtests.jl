using Symplectic
using Test
using LinearAlgebra

const files = (
    "symp",
    "symp_generators",
    # "color_sets",
    # "interference_based_sequence",
    # "teleportation_based_generator",
    # "cayley_transform"
)

@testset "Symplectic.jl" begin    
    @testset "$(titlecase(f))" for f in files
        include("$f.jl")
    end
end

# teleportationBasedSymplecticControl( Symplectic.randomGenericSymp(2), [1, 2], [1, 2]  )
# randomizeSymp( beamSplitter(1), 2  )