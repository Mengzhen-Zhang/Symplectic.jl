using Test
using Symplectic

function test_adaptive(S, in, out)
    F = feedforward(S, in, out)
    A = adaptiveMeasurement(F, out, 2)
    out1 = vcat([[2*o-1, 2*o] for o in out]...)
    in1 = vcat([[2*i-1, 2*i] for i in in]...)
    (A * S)[out1, in1] ≈ teleport(S, in, out)
end

@testset "lib" begin
    @test test_adaptive(bs(π/3), [1], [2])
    @test test_adaptive(bs(π/6), [1], [2])
    @test test_adaptive(bs(π/6), [1], [1])
    @test test_adaptive(bs(π/4), [2], [2])
end