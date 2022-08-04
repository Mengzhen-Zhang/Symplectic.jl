export localSympFromQuad, getDecoupleSequence, getInterferenceBasedSequence, getSequence

# From (source quadrautre, target quadratre) to Local Symplectic Matrix transforming source to Ω target
↑(m, i) = 2*m - (i % 2)
↓(m, i) = 2*m -1 + (i  % 2)
getPattern(v::Vector)=(n->n>tolerance).( (diag∘sqrt∘Diagonal∘(x->x'*x)∘(v->reshape(v, 2, :)))(v) )
function Lkmij(u, v, m, i, j)
    if norm(v[2*m-1:2*m])<tolerance
        return i==j ? 1 : 0
    else
        return (-1)^(j+1)*v[↑(m, i)]*u[↓(m, j)]/( v[2*m]^2 + v[2*m-1]^2 ) + (-1)^i*v[↓(m, i)]*u[↑(m, j)]/( u[2*m]^2 + u[2*m-1]^2 )
    end
end
function localSympFromQuad(u::Vector, v::Vector; lambda::Real = 1)::Symp
    if length(u) % 2 !=0 || length(u) != length(v) || getPattern(u) != getPattern(v) 
        return error("input not valid") 
    end
    n = length(u) ÷ 2
    if n == 1
        return Symp([ Lkmij(lambda * u, v, 1, i, j) for i in 1:2, j in 1:2 ])
    else
        return ⊕([ Symp([ Lkmij(lambda * u, v, m, i, j) for i in 1:2, j in 1:2 ]) for m in 1:n ]...)
    end
end
# Get Local Operations L3 L2 L1 for the Decoupling Sequence S4 L3 S3 L2 S2 L1 S1
function getDecoupleSequence(S4::Symp, S3::Symp, S2::Symp, S1::Symp, m::Int=1)::Vector{Symp}
    # if !(all(x->isGeneric(x), [S1, S2, S3, S4]) ) return error("not generic") end
    L1 = localSympFromQuad(S1.S[:,2*m-1], S2.S[2*m-1,:])
    L3 = localSympFromQuad(S3.S[:,2*m-1], S4.S[2*m-1,:])
    T1 = S2 * L1 * S1
    T2 = S4 * L3 * S3
    u = T1.S[:,2*m]
    v = T2.S[2*m,:] + sum( (k->T1.S[k,2*m-1]*T1.S[k,2*m]-T2.S[2*m-1,k]*T2.S[2*m,k]).(1:length(u))  ) * T2.S[2*m-1, :]
    L2 = localSympFromQuad(u, v; lambda = lambda)
    return [S4, L3, S3, L2, S2, L1, S1]
end

# Get local Operations L16, L15, ..., L1 for the interfernce-basded sequence LR S L15 S ... S L1 S, with a given 4×4 target Symplectic Matrix S⊙
function getSequence(Ss::Vector{Symp{T}}, modeToDecouple::Int)::Vector{Symp} where T
    if length(Ss) == 4
        return getDecoupleSequence(Ss..., modeToDecouple)
    else
        step = length(Ss) ÷ 4
        Rs = [ getSequence(Ss[i:i+step-1], modeToDecouple-1 ) for i in 1:step:length(Ss) ]

        Ls = getDecoupleSequence([ *(RLst...) for RLst in Rs]..., modeToDecouple)
        return vcat([ (i%2==0) ? [ Ls[i] ] : Rs[(i+1) ÷ 2] for i in 1:7 ]...)
    end
end

function getInterferenceBasedSequence(S::Symp, ST::Symp)::Vector{Symp}
    if !(isGeneric(S))  return error("not generic") end
    Ss = [S for i in 1:4^nModes(ST)]
    Ss[end] =  Ss[end] * (ST ⊕ Id(nModes(S) - nModes(ST)))^-1
    Ss = getSequence(Ss, nModes(ST))
    LR = Symp( ( *(Ss...).S[1:2*nModes(ST), 1:2*nModes(ST)]) )^-1  ⊕ Id(nModes(S) - nModes(ST)) 
    Ss[end] =  Ss[end] * (ST ⊕ Id(nModes(S) - nModes(ST)))
    return vcat([ LR ] , Ss)
end



# Interference-Based Sequence with Automatic Local Randomization
function getSequenceWithRandomization(Ss::Vector{Symp{T}}, modeToDecouple::Int)::Vector{Symp} where T
    if length(Ss) == 4
        return getDecoupleSequence(Ss..., modeToDecouple)
    else
        step = length(Ss) ÷ 4
        Rs = [ getSequence(Ss[i:i+step-1], modeToDecouple-1 ) for i in 1:step:length(Ss) ]

        Ls = getDecoupleSequence([ *(RLst...) for RLst in Rs]..., modeToDecouple)
        return vcat([ (i%2==0) ? [ Ls[i] ] : Rs[(i+1) ÷ 2] for i in 1:7 ]...)
    end
end
function getIQSWithRandomization(S::Symp, ST::Symp)::Vector{Symp}
    if !(isGeneric(S))  return error("not generic") end
    Ss = [S for i in 1:4^nModes(ST)]
    Ss[end] =  Ss[end] * (ST ⊕ Id(nModes(S) - nModes(ST)))^-1
    Ss = getSequenceWithRandomization(Ss, nModes(ST))
    LR = Symp( ( *(Ss...).S[1:2*nModes(ST), 1:2*nModes(ST)]) )^-1  ⊕ Id(nModes(S) - nModes(ST)) 
    Ss[end] =  Ss[end] * (ST ⊕ Id(nModes(S) - nModes(ST)))
    return vcat([ LR ] , Ss)
end