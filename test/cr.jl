using Symplectic
using LinearAlgebra
using Optim

cr = CoupledResonators()

cr = addPassive(cr, 0.0, 1, 1)
cr = addPassive(cr, 1.0, 1, 1)

cr = addPassive(cr, 1.0, 1, 2 )
cr = addActive(cr, 1.0, 1, 2 )

S = scatteringMatrix( 0.0, cr  )

r = [1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 1 0 0;
    0.0 0 0 0 0 0;
    0 0.0 0 0 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1]

l = [1 0 0 0 0.0 0 0 0;
     0 1 0 0 0 0.0 0 0;
     0 0 1 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
     0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 1]


nonSymplecticity(l * S * r)
l * Ω * r