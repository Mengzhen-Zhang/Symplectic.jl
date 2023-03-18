using Symplectic
using LinearAlgebra
# using Optim
# using Plots
# using LaTeXStrings
# using Plots.PlotMeasures

cr = CoupledResonators()

cr = addPassive(cr, 1.0, 1, 1)
cr = addPassive(cr, 1.0, 2, 2)
cr = addPassive(cr, 1.0, 3, 3)
cr = addPassive(cr, 1.0, 4, 4)
cr = addPassive(cr, 1.0, 5, 5)

cr = addPassive(cr, 1.0, 1, 2 )
cr = addActive(cr, 1.0, 1, 2 )

cr = addPassive(cr, 1.0, 1, 3 )
cr = addActive(cr, 1.0, 1, 3 )

cr = addPassive(cr, 1.0, 1, 4 )
cr = addActive(cr, 1.0, 1, 4 )

cr = addPassive(cr, 1.0, 5, 2 )
cr = addActive(cr, 1.0, 5, 2 )

cr = addPassive(cr, 1.0, 5, 3 )
cr = addActive(cr, 1.0, 5, 3 )

cr = addPassive(cr, 1.0, 5, 4 )
cr = addActive(cr, 1.0, 5, 4 )

cr = add

S = scatteringMatrix( 1.0, cr  )

# const bs = beamSplitter

# function evaluation(θs::Vector)
#     S = solution(θs)
#     T = bs(π / 2) # [0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0] # [0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0] # bs(π/2) # I(4) # bs(π/3) * Diagonal([2, 1/2, 1/3, 3]) * bs(π/3)
#     return sqrt(tr((S - T) * transpose(S - T)))
# end

# function solution(θs::Vector)
#     return teleportation(phaseShifting([1, 1, θs[1:8]...]...)*Matrix(S)*phaseShifting([1, 1, θs[9:16]...]...), [1, 2], [1, 2])
# end

# res = optimize(evaluation,
#                 2*π*rand(16),
#                 method = LBFGS(), 
#                 iterations = 1000000,
#                 autodiff = :forward)
# round.(solution(Optim.minimizer(res)), digits = 4)
# Optim.minimizer(res)

# function transmission(θ)
#     res = teleportation(
#         phaseShifting(
#             [1, 1, θ, θ]...)*Matrix(S)*phaseShifting([1, 1, θ + π/3, θ + π/6]...), 
#             [1, 2],
#             [1, 2])
#     return det(res[3:4, 1:2])
# end

# plot(transmission, -π, π, 
#         ylims = (-2, 2),
#         xguide = L"$\theta$",
#         xguide_position = :center,
#         xticks = (π*[-1, -2/3,  -1/3,  0,  1/3,  2/3, 1], 
#                     [L"\pi", L"-\frac{2\pi}{3}",  L"-\frac{\pi}{3}",  L"0",  L"\frac{\pi}{3}", L"\frac{2\pi}{3}", L"\pi"]),
#         xtickfontsize = 11,
#         yticks = ([-2, -1, 0, 1, 2],
#                  [L"-2", L"-1", L"0", L"1", L"2"]),
#         ytickfontsize = 11,
#         label = "Transmission",
#         bottom_margin = 10px )

# savefig("~/test.pdf")