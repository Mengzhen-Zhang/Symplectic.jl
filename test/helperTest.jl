using Symplectic
using Test

@testset "CustomFunction" begin
    func1 = (x) -> x + 1
    func2 = (x, y) -> x * y
    f1 = CustomFunction(1, func1)
    f2 = CustomFunction(2, func2)

    @test begin
        f1(1) == func1(1) && f2(3, 4) == func2(3, 4)
    end

    @test begin
        (f1 + f2)(1, 3, 4) == func1(1) + func2(3, 4)
    end

    @test begin
        (f1 * f2)(1, 3, 4) == func1(1) * func2(3, 4)
    end
end