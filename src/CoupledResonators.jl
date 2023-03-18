using LinearAlgebra:
        Diagonal, copyto!, diag, I
using SkewLinearAlgebra

export CoupledResonators
export addModes, addActive, addPassive
export addGammaEx, addGammaIn
export scatteringMatrix
export dilatedScatteringMatrix

struct CoupledResonators
    yr::AbstractMatrix
    yi::AbstractMatrix
    wr::AbstractMatrix
    wi::AbstractMatrix
    γin::AbstractMatrix
    γex::AbstractMatrix
end

function CoupledResonators()
    return CoupledResonators(
        zeros(0, 0),
        zeros(0, 0),
        zeros(0, 0),
        zeros(0, 0),
        Diagonal{Real}([]),
        Diagonal{Real}([]),
     )
end

modes(cr::CoupledResonators) = size(cr.yr, 1)

function addModes(cr::CoupledResonators, l::Int)
    if l == 0
        return cr
    end
    return CoupledResonators(
    Matrix(cr.yr ⊕ zeros(l, l)),
    Matrix(cr.yi ⊕ zeros(l, l)),
    Matrix(cr.wr ⊕ zeros(l, l)),
    Matrix(cr.wi ⊕ zeros(l, l)),
    Diagonal( [ diag(cr.γin)..., zeros(l)... ] ),
    Diagonal( [ diag(cr.γex)..., ones(l)... ] ) 
    )
end

# add g a⁺ⱼaₖ + conj(g) aⱼa⁺ₖ
function addPassive(cr::CoupledResonators, 
    g::Number, j::Int, k::Int)
    m = size(cr.yr, 1)
    cr = addModes(cr, max(j - m, k - m, 0 ))
    yr = cr.yr
    yi = cr.yi
    wi = cr.wi
    wr = cr.wr
    γin = cr.γin
    γex = cr.γex
    yr[j, k] = real(g)
    yr[k, j] = real(g)
    yi[j, k] = imag(g)
    yi[k, j] = -imag(g)
    return CoupledResonators(yr, yi, wr, wi, γin, γex)
end

# add g a⁺ⱼa⁺ₖ + conj(g) aⱼaₖ
function addActive(cr::CoupledResonators, 
    g::Number, j::Int, k::Int)
    m = size(cr.yr, 1)
    cr = addModes(cr, max(j - m, k - m, 0 ))
    yr = cr.yr
    yi = cr.yi
    wi = cr.wi
    wr = cr.wr
    γin = cr.γin
    γex = cr.γex
    wr[j, k] = real(g) / 2
    wr[k, j] = real(g) / 2
    wi[j, k] = imag(g) / 2
    wi[k, j] = imag(g) / 2
    return CoupledResonators(yr, yi, wr, wi, γin, γex)
end

function addGammaIn(cr::CoupledResonators,
    γ::Number, l::Int)
    m = size(cr.yr, 1)
    cr = addModes(cr, max(l - m, 0 ))
    yr = cr.yr
    yi = cr.yi
    wi = cr.wi
    wr = cr.wr
    γin = cr.γin
    γex = cr.γex
    γin[l, l] = γ
    return CoupledResonators(yr, yi, wr, wi, γin, γex)
end

function addGammaEx(cr::CoupledResonators,
    γ::Number, l::Int)
    m = size(cr.yr, 1)
    cr = addModes(cr, max(l - m, 0 ))
    yr = cr.yr
    yi = cr.yi
    wi = cr.wi
    wr = cr.wr
    γin = cr.γin
    γex = cr.γex
    γex[l, l] = γ
    return CoupledResonators(yr, yi, wr, wi, γin, γex)
end

const e0 = I(4)

const e1 = [0 0 0 1;
            0 0 1 0;
            0 1 0 0;
            1 0 0 0]

const e2 = [0 0 1 0;
            0 0 0 -1;
            1 0 0 0;
            0 -1 0 0]

const e3 = [1 0 0 0;
            0 1 0 0;
            0 0 -1 0;
            0 0 0 -1]

function scatteringMatrix(ω, yr, yi, wr, wi, γin, γex)
    m = size(yr, 1)
    Θ = ω * inv(γex + γin)
    Yr = inv(sqrt.(γex)) * yr * inv(sqrt.(γex))
    Yi = inv(sqrt.(γex)) * yi * inv(sqrt.(γex))
    Wr = inv(sqrt.(γex)) * wr * inv(sqrt.(γex))
    Wi = inv(sqrt.(γex)) * wi * inv(sqrt.(γex))
    Γ = inv(sqrt.(γex)) * (γex + γin) * inv(sqrt.(γex)) / 2
    S = e0 ⊗ I(m) - inv(
          e0 ⊗ (Yi + Γ) -
          e1 ⊗ Wr +
          e2 ⊗ Wi -
          (e1 * e2) ⊗ Yr +
          (e1 * e2 * e3) ⊗ Θ
    )
    order = vcat([ [i, i + m] for i in 1:m]...,  [ [i + 2*m, i + 3*m] for i in 1:m]...)
    return S[order, order]
end
function scatteringMatrix(ω::Number, cr::CoupledResonators)
    yr = cr.yr
    yi = cr.yi
    wr = cr.wr
    wi = cr.wi
    γin = cr.γin
    γex = cr.γex
    return scatteringMatrix(ω, yr, yi, wr, wi, γin, γex)
end

function dilatedScatteringMatrix(ω::Number, cr::CoupledResonators)
    S = scatteringMatrix(ω, cr)
    return dilate(S)
end