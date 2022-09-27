export teleportation

"""
    teleportation(S::AbstractMatrix, inModes::Vector, outModes::Vector)

Caculate Eq. (4) in *PHYSICAL REVIEW LETTERS 120, 020502 (2018)*, with modes not contained in
`inModes` squeezed along their Q-quadratures and modes not contained in `outModes` measured
along their P-quadratures.
"""
function teleportation(S::AbstractMatrix, inModes::Vector, outModes::Vector)
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
       SOutIn = S[Out, In]
       SHmIn = S[Hm, In]
       SOutUsq = S[Out, Usq]
       SHmUsq = S[Hm, Usq]
       return SOutIn - SOutUsq * Matrix(SHmUsq)^-1 * SHmIn
end
