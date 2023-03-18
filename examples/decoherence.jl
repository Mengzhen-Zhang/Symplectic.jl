using Symplectic, LinearAlgebra

begin
    cr = CoupledResonators()
    cr = addPassive(cr, 1.0, 1, 1)
    cr = addPassive(cr, 1.0, 2, 2)
    cr = addPassive(cr, 1.0, 1, 2 )
    cr = addActive(cr, 1.0, 1, 2 )
    cr = addGammaIn(cr, 0.1, 1)    
end

Sd = dilatedScatteringMatrix(1.0, cr)


