tolerance = 10^-6

S = randomGenericSymp(4)

@test size(S) == (8, 8)
@test typeof(S[1, 1]) == BigFloat
@test typeof(S) <: Symp
@test typeof(S) <: AbstractMatrix
@test typeof(Symp(S)) <: Symp
@test typeof(Matrix(S)) <: AbstractMatrix
@test isGeneric(S) == true
@test nModes(S) == 4
@test typeof(transpose(S)) <: Symp
@test typeof(S') <: Symp
@test typeof(-S) <: Symp
@test typeof(+S) <: Symp
@test isapprox( inv(S) * S, Id(4); atol = tolerance )
@test Omega(2) == Symp([0 1 0 0; -1 0 0 0; 0 0 0 1; 0 0 -1 0])
@test Id(2) == Symp([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1])
@test dsum(Id(2), Id(2), Id(2)) == Id(6)
@test S^-1 == inv(S)