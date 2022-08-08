using Symplectic
import Optim
import LinearAlgebra: det


function target(tel)
    s = tel[1:2, 3:4]
    return (det(s))^2
end

res = findSolutionFromFile("./test/test_optimize.json", target)

sc = buildCircuitFromFile("./test/test_optimize.json")

so = reduce(*, sc.circuit)
fin = teleportation(so(Optim.minimizer(res)...), [1, 2], [1, 2])

det(fin[3:4, 1:2]) < 1e-9
det(fin[1:2, 3:4]) < 1e-9
Symplectic.nonSymplecticity(fin) < 1e-9
