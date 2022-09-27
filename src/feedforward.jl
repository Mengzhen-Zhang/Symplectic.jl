export feedforward

"""
    feedforward(S::AbstractMatrix, inModes::Vector, outModes::Vector)

Caculate Eq. (3) in *PHYSICAL REVIEW LETTERS 120, 020502 (2018)*, with modes not contained in
`inModes` squeezed along their Q-quadratures and modes not contained in `outModes` measured
along their P-quadratures.
"""
function feedforward(S::AbstractMatrix, inModes::Vector, outModes::Vector)
       if length(inModes) == 0
              return S
       end
       allModes = [1:(size(S, 1)÷2);]
       In = vcat([[2 * i - 1, 2 * i] for i in inModes]...)
       Out = vcat([[2 * i - 1, 2 * i] for i in outModes]...)
       ancModes = [i for i in allModes if !(i in inModes)]
       idlModes = [i for i in allModes if !(i in outModes)]
       Usq = 2 * ancModes
       Hm = 2 * idlModes
       SOutUsq = S[Out, Usq]
       SHmUsq = S[Hm, Usq]
       return - SOutUsq * Matrix(SHmUsq)^-1
end
