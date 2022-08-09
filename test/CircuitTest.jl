using Symplectic
using Test

import Optim
import LinearAlgebra: det

@testset "regularize" begin
    regularize = Symplectic.regularize

    @test regularize(3, PhaseShifting, Any["_", 1.0, 2.0]) == Any["_", 1.0, 2.0]
    @test regularize(3, BeamSplitter, Any["_", 1, 2]) == Any["_", 1, 2, 3]
    @test regularize(1, Custom, Any["_", "_", "_", "_"]) == Any["_", "_", "_", "_", 2]
end

@testset "buildCircuit" begin
    buildCircuitFromFile = Symplectic.buildCircuitFromFile

    sc = buildCircuitFromFile("./test/test.json")
    @test length(sc.circuit) == 9
    @test length(sc.inModes) == 0
    @test length(sc.outModes) == 0

    ps = phaseShifting(1.0, 2.0, 3.0, 4.0)
    @test sc.circuit[9](1.0, 2.0, 3.0, 4.0) == ps
    @test sc.circuit[8](1.0) == ps
    @test sc.circuit[7]() == beamSplitter(1.0, 1, 2, 4)
    @test sc.circuit[6]() == beamSplitter(2.0, 2, 3, 4)
    @test sc.circuit[5]() == beamSplitter(3.0, 3, 4, 4)
    @test sc.circuit[4](1.0) == beamSplitter(1.0, 1, 4, 4)
    @test sc.circuit[3](repeat([1.0], 64)...) == reshape(repeat([1.0], 64), (8, 8))
    @test sc.circuit[2]() == reshape(
        [
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0
        ],
        (8, 8)
    )
    @test sc.circuit[1]() == phaseShifting(2.0, 3.0, 4.0, 1.0)
end

@testset "buildConstraint" begin
    sc = Symplectic.buildCircuitFromFile("./test/test.json")
    @test Symplectic.buildConstraint(sc).l == 70
end
