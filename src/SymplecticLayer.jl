using ChainRulesCore

export SymplecticLayer, SymplecticCircuit, symplecticLayer
export GaussianChannelLayer, GaussianChannelCircuit, gaussianChannelLayer

struct SymplecticLayer
    num_of_args::Integer
    num_of_modes::Integer
    is_symplectic::Bool
    S::Function
end

(sl::SymplecticLayer)(x::AbstractVector) = sl.S(x)

Base.:*(sl::SymplecticLayer) = sl
Base.:*(sl1::SymplecticLayer, sl2::SymplecticLayer) = SymplecticLayer(
        sl1.num_of_args + sl2.num_of_args,
        sl1.num_of_modes,
        sl1.is_symplectic && sl2.is_symplectic,
        x -> sl1.S(x[begin : sl1.num_of_args]) * sl2.S(x[sl1.num_of_args + 1 : end])
    )

struct SymplecticCircuit
    layers::AbstractVector{SymplecticLayer}
    is_adaptive::Bool
    inModes::AbstractVector{Integer}
    outModes::AbstractVector{Integer}
    ancModes::AbstractVector{Integer}
    idlModes::AbstractVector{Integer}
    envModes::AbstractVector{Integer}
    num_of_args::Integer
    num_of_modes::Integer
    S::Function
    function SymplecticCircuit(layers, is_adaptive, inModes, 
        outModes, ancModes, idlModes, num_of_modes)
        envModes = [m for m in 1:num_of_modes if m ∉ ancModes && m ∉ inModes]
        preps = SymplecticLayer(
            length(ancModes),
            num_of_modes,
            true,
            x -> phaseShifting([m ∈ ancModes ? x[findfirst(==(m), ancModes)] : 0 for m in 1:num_of_modes]...)
        )
        postps = SymplecticLayer(
            length(idlModes),
            num_of_modes,
            true,
            x -> phaseShifting([m ∈ idlModes ? x[findfirst(==(m), idlModes)] : 0 for m in 1:num_of_modes]...)
        )
        len_of_F = 2 * length(outModes) * length(idlModes)
        raw = *(layers...)
        A = SymplecticLayer(len_of_F,
            num_of_modes,
            true,
            x -> adaptiveMeasurement(
                [reshape(x[begin:len_of_F], (2*length(outModes), length(idlModes))); 
                 zeros(2*length(envModes), length(idlModes))],
                vcat(outModes, envModes), num_of_modes)
        )
        final = if is_adaptive
            A * postps * raw * preps
        else
            raw
        end
        num_of_args, S = final.num_of_args, final.S
        return new(layers, is_adaptive, inModes,
        outModes, ancModes, idlModes, envModes,
        num_of_args, num_of_modes, S)
    end
end

(sc::SymplecticCircuit)(x) = sc.S(x)

function symplecticLayer(ω::Float64, cr::CoupledResonators)
    return SymplecticLayer(
        0,
        size(cr.γex, 1) * 2,
        true,
        x -> dilatedScatteringMatrix(ω, cr)
    )
end

struct GaussianChannelLayer 
    num_of_args::Integer
    num_of_modes::Integer
    T::Function
    N::Function
end

struct GaussianChannelCircuit
    layers::AbstractVector{GaussianChannelLayer}
    is_adaptive::Bool
    inModes::AbstractVector{Integer}
    outModes::AbstractVector{Integer}
    ancModes::AbstractVector{Integer}
    idlModes::AbstractVector{Integer}
    num_of_args::Integer
    num_of_modes::Integer
    squeezing::Float64
    T::Function
    N::Function
    function GaussianChannelCircuit(layers, is_adaptive, inModes, outModes,
        ancModes, idlModes, squeezing)
        num_of_modes = length(inModes) + length(ancModes)
        preps = SymplecticLayer(
            length(ancModes),
            num_of_modes,
            true,
            x -> phaseShifting([m ∈ ancModes ? x[findfirst(==(m), ancModes)] : 0 for m in 1:num_of_modes]...)
        )
        postps = SymplecticLayer(
            length(idlModes),
            num_of_modes,
            true,
            x -> phaseShifting([m ∈ idlModes ? x[findfirst(==(m), idlModes)] : 0 for m in 1:num_of_modes]...)
        )
        len_of_F = 2 * length(outModes) * length(idlModes)
        A = SymplecticLayer(len_of_F,
            num_of_modes,
            true,
            x -> adaptiveMeasurement(
                reshape(x[begin:len_of_F], (2*length(outModes), length(idlModes))),
                outModes, num_of_modes)
        )
        raw = *(layers...)
        final = if is_adaptive
            gl = A * postps * raw * preps
            in_indices = vcat([[2*m-1, 2*m] for m in inModes]...)
            out_indices = vcat([[2*m-1, 2*m] for m in outModes]...)
            GaussianChannelLayer(gl.num_of_args, length(outModes),
                x -> begin
                    T = gl.T(x[begin:length(x)])
                    T[out_indices, in_indices]
                end,
                x -> begin
                    T = gl.T(x[begin:length(x)])
                    N = gl.N(x[begin:length(x)])
                    anc_indices = vcat([[2*m-1, 2*m] for m in ancModes]...)
                    (T[:, anc_indices]*squeezedVacuum(fill(squeezing, length(ancModes)))*transpose(T[:, anc_indices]) 
                    + N)[out_indices, out_indices]
                end)
        else
            gl = raw
            in_indices = vcat([[2*m-1, 2*m] for m in inModes]...)
            out_indices = vcat([[2*m-1, 2*m] for m in outModes]...)
            GaussianChannelLayer(gl.num_of_args, gl.num_of_modes,
                x -> gl.T(x[begin:length(x)])[out_indices, in_indices],
                x -> begin
                    T = gl.T(x[begin:length(x)])
                    N = gl.N(x[begin:length(x)])
                    anc_indices = vcat([[2*m-1, 2*m] for m in ancModes]...)
                    (T[:, anc_indices]*transpose(T[:, anc_indices]) + N)[out_indices, out_indices]
                end)
        end
        num_of_args, T, N = final.num_of_args, final.T, final.N
        return new(layers, is_adaptive, inModes, outModes, ancModes, idlModes,
        num_of_args, num_of_modes, squeezing, T, N)
    end
end

(gl::GaussianChannelLayer)(x::AbstractVector) = (gl.T(x), gl.N(x))

Base.:*(gl::GaussianChannelLayer) = gl
Base.:*(gl1::GaussianChannelLayer, gl2::GaussianChannelLayer) = GaussianChannelLayer(
            gl1.num_of_args + gl2.num_of_args,
            gl1.num_of_modes,
            x -> begin
                T1 = gl1.T(x[begin : gl1.num_of_args])
                T2 = gl2.T(x[gl1.num_of_args + 1 : end])
                T1 * T2
            end,
            x -> begin
                N1 = gl1.N(x[begin : gl1.num_of_args])
                N2 = gl2.N(x[gl1.num_of_args + 1 : end])
                T2 = gl2.T(x[gl1.num_of_args + 1 : end])
                T1 * N2 * transpose(T1) + N1
            end
        )
Base.:*(sl::SymplecticLayer, gl::GaussianChannelLayer) = GaussianChannelLayer(
    sl.num_of_args + gl.num_of_args,
    sl.num_of_modes,
    x -> begin
        S = sl.S(x[begin : sl.num_of_args])
        T = gl.T(x[sl.num_of_args + 1 : end])
        S * T
    end,
    x -> begin
        S = sl.S(x[begin : sl.num_of_args])
        N = gl.N(x[sl.num_of_args + 1 : end])
        S * N * transpose(S)
    end
)

Base.:*(gl::GaussianChannelLayer, sl::SymplecticLayer) = GaussianChannelLayer(
    sl.num_of_args + gl.num_of_args,
    gl.num_of_modes,
    x -> begin
        T = gl.T(x[begin : gl.num_of_args])
        S = sl.S(x[gl.num_of_args + 1 : end])
        T * S
    end,
    x -> begin
        N = gl.N(x[begin : gl.num_of_args])
        N
    end
)
 
function gaussianChannelLayer(ω::Float64, cr::CoupledResonators)
    num_of_modes = size(cr.γex, 1) * 2
    Sd = dilatedScatteringMatrix(ω, cr)
    return GaussianChannelLayer(
        0, 
        num_of_modes,
        x -> Sd[begin:2*num_of_modes, begin:2*num_of_modes],
        x -> Sd[begin:2*num_of_modes, 2*num_of_modes+1:end] * transpose(Sd[begin:2*num_of_modes, 2*num_of_modes+1:end])
    )
end