using LinearAlgebra: I

@testset "toSympType" begin
    @test toSympType("beam_splitter") == BeamSplitter
    @test toSympType("phase_shifting") == PhaseShifting
    @test toSympType("custom") == Custom
    @test toSympType("random_name") == Undef
end


@testset "SymplecticOperation" begin
    bs1 = SymplecticOperation(3, BeamSplitter, Any[1.0, 1, 2, 3])
    bs2 = SymplecticOperation(3, BeamSplitter, Any["_", 1, 2, 3])
    ps1 = SymplecticOperation(4, PhaseShifting, Any[1.0, 2.0, 3.0, 4.0])
    ps2 = SymplecticOperation(4, PhaseShifting, Any[1.0, "_", "_", 4.0])
    ps3 = SymplecticOperation(4, PhaseShifting, Any["_", "_", "_", "_"])
    custom1 = SymplecticOperation(1, Custom, Any[1.0, 2.0, 3.0, 4.0, 2])
    custom2 = SymplecticOperation(1, Custom, Any[1.0, "_", "_", 4.0, 2])
    custom3 = SymplecticOperation(1, Custom, Any["_", "_", "_", "_", 2])

    @test begin
        bs1() == beamSplitter(1.0, 1, 2, 3)
    end

    @test begin
        bs2(1.0) == beamSplitter(1.0, 1, 2, 3)
    end

    @test begin
        ps1(1.0) == phaseShifting(1.0, 2.0, 3.0, 4.0)
    end

    @test begin
        ps2(2.0, 3.0) == phaseShifting(1.0, 2.0, 3.0, 4.0)
    end

    @test begin
        ps3(1.0, 2.0, 3.0, 4.0) == phaseShifting(1.0, 2.0, 3.0, 4.0)
    end

    @test begin
        custom1() == reshape([1.0, 2.0, 3.0, 4.0], (2, 2))
    end

    @test begin
        custom2(2.0, 3.0) == reshape([1.0, 2.0, 3.0, 4.0], (2, 2))
    end

    @test begin
        custom3(1.0, 2.0, 3.0, 4.0) == reshape([1.0, 2.0, 3.0, 4.0], (2, 2))
    end

    s1 = SymplecticOperation(2, BeamSplitter, Any["_", 1, 2, 2])
    s2 = SymplecticOperation(3, BeamSplitter, Any["_", 2, 3, 3])

    @test begin
        (s2 ⊕ s1)(1.0, 1.0) == beamSplitter(1.0, 2, 3, 3) ⊕ beamSplitter(1.0, 1, 2, 2)
    end

    @test begin
        n1 = 2
        (s2 ⊕ n1)(1.0, 1.0) == beamSplitter(1.0, 2, 3, 3) ⊕ I(4)
    end

    @test begin
        n2 = 3
        (n2 ⊕ s1)(1.0, 1.0) == I(6) ⊕ beamSplitter(1.0, 1, 2, 2)
    end

    @test begin
        (s2 * s1)(1.0, 1.0) == beamSplitter(1.0, 2, 3, 3) * beamSplitter(1.0, 1, 2, 3)
        (s1 * s2)(1.0, 1.0) == beamSplitter(1.0, 1, 2, 3) * beamSplitter(1.0, 2, 3, 3)
    end

    @test begin
        inv(s2)(1.0) == -Ω * transpose(beamSplitter(1.0, 2, 3, 3)) * Ω
    end

    @test begin
        teleportation(s2, [1, 2], [1, 2])(1.0) == teleportation(beamSplitter(1.0, 2, 3, 3), [1, 2], [1, 2])
    end

    @test begin
        Symplectic.nonSymplecticity(custom2)(2.0, 3.0) == Symplectic.nonSymplecticity(reshape([1.0, 2.0, 3.0, 4.0], (2, 2)))
    end

end