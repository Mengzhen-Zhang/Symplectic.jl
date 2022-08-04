export  teleportationBasedSympGenerator

# Calculate the effective Symplectic Matrix output of the adaptive Gaussian control scheme
function teleportationBasedSympGenerator(S::Symp, in_modes::Vector, out_modes::Vector)::Symp
    all_modes = [1:nModes(S);]
    in_ = vcat( [ [2*i - 1, 2 * i] for i in in_modes ]... )
    out = vcat( [ [2*i - 1, 2 * i] for i in out_modes ]... )
    ancilliary_modes = [i for i in all_modes if !(i in in_modes)]
    idle_modes = [i for i in all_modes if !(i in out_modes)]
    usq = 2*ancilliary_modes
    hm = 2*idle_modes
    S_out_in = S.S[ out, in_  ]
    S_hm_in = S.S[ hm, in_ ]
    S_out_usq = S.S[ out, usq ]
    S_hm_usq = S.S[ hm, usq ]
    return S_out_in - S_out_usq * S_hm_usq^-1 * S_hm_in
end