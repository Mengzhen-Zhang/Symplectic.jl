using Symplectic
using Test

@testset "toSympType" begin
    @test toSympType("beam_splitter") == BeamSplitter
    @test toSympType("phase_shifting") == PhaseShifting
    @test toSympType("custom") == Custom
    @test toSympType("random_name") == Undef
end

@testset "funcBuild" begin
    @test begin
        f = funcBuild(beamSplitter, Any[1.0, 1, 2, 3])
        f() == beamSplitter(1.0, 1, 2, 3)
    end

    @test begin
        f = funcBuild(beamSplitter, Any["_", 1, 2, 3])
        f(1.0) == beamSplitter(1.0, 1, 2, 3)
    end

    @test begin
        f = funcBuild(phaseShifting, Any[1.0, 2.0, 3.0, 4.0])
        f() == phaseShifting(1.0, 2.0, 3.0, 4.0)
    end

    @test begin
        f = funcBuild(phaseShifting, Any["_", "_", "_", "_"])
        f(1.0, 2.0, 3.0, 4.0) == phaseShifting(1.0, 2.0, 3.0, 4.0)
    end
    
    @test begin
        f = funcBuild(phaseShifting, Any[1.0, "_", "_", 4.0])
        f(2.0, 3.0) == phaseShifting(1.0, 2.0, 3.0, 4.0)
    end

    @test begin
        f = funcBuild(customReshape, Any[1.0, 2.0, 3.0, 4.0, 2])
        f() == reshape([1.0, 2.0, 3.0, 4.0], (2, 2))
    end

    @test begin
        f = funcBuild(customReshape, Any[1.0, "_", 3.0, 4.0, 2])
        f(2.0) == reshape([1.0, 2.0, 3.0, 4.0], (2, 2))
    end

    @test begin
        f = funcBuild(customReshape, Any["_", "_", "_", "_", 2])
        f(1.0, 2.0, 3.0, 4.0) == reshape([1.0, 2.0, 3.0, 4.0], (2, 2))
    end
end

@testset "SymplecticOperation" begin
    @test begin
        so = SymplecticOperation(3, "beam_splitter", Any[1.0, 1, 2])
        so.type == BeamSplitter && so.Op() == beamSplitter(1.0, 1, 2, 3)
    end

    @test begin
        so = SymplecticOperation(3, "beam_splitter", Any["_", 1, 2])
        so.type == BeamSplitter && so.Op(1.0) == beamSplitter(1.0, 1, 2, 3)
    end

    @test begin
        so = SymplecticOperation(4, "phase_shifting", Any[1.0, 2.0, 3.0, 4.0])
        so.type == PhaseShifting && so.Op(1.0) == phaseShifting(1.0, 2.0, 3.0, 4.0)
    end

    @test begin
        so = SymplecticOperation(4, "phase_shifting", Any[1.0, "_", "_", 4.0])
        so.type == PhaseShifting && so.Op(2.0, 3.0) == phaseShifting(1.0, 2.0, 3.0, 4.0)
    end

    @test begin
        so = SymplecticOperation(4, "phase_shifting", Any["_"])
        so.type == PhaseShifting && so.Op(1.0, 2.0, 3.0, 4.0) == phaseShifting(1.0, 2.0, 3.0, 4.0)
    end

    @test begin
        so = SymplecticOperation(1, "custom", Any[1.0, 2.0, 3.0, 4.0])
        so.type == Custom && so.Op() == reshape([1.0, 2.0, 3.0, 4.0], (2, 2))
    end

    @test begin
        so = SymplecticOperation(1, "custom", Any[1.0, "_", "_", 4.0])
        so.type == Custom && so.Op(2.0, 3.0) == reshape([1.0, 2.0, 3.0, 4.0], (2, 2))
    end

    @test begin
        so = SymplecticOperation(1, "custom", Any["_"])
        so.type == Custom && so.Op(1.0, 2.0, 3.0, 4.0) == reshape([1.0, 2.0, 3.0, 4.0], (2, 2))
    end

    @test begin
        s1 = SymplecticOperation(3, "beam_splitter", Any["_", 1, 2])
        s2 = SymplecticOperation(3, "beam_splitter", Any["_", 2, 3])
        (s2 * s1).Op(1.0,1.0) == beamSplitter(1.0, 2, 3, 3) * beamSplitter(1.0, 1, 2, 3) 
    end

end