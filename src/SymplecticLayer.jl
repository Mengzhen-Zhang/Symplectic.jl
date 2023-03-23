using ChainRulesCore

export SymplecticLayer, SymplecticCircuit, symplecticLayer

struct SymplecticLayer
    num_of_args::Integer
    f::Function
end

(sl::SymplecticLayer)(x::AbstractVector) = sl.f(x)

Base.:*(sl::SymplecticLayer) = sl
Base.:*(sl1::SymplecticLayer, sl2::SymplecticLayer) = SymplecticLayer(
        sl1.num_of_args + sl2.num_of_args,
        x -> sl1.f(x[begin : sl1.num_of_args]) * sl2.f(x[sl1.num_of_args + 1 : end])
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
    f::Function
    function SymplecticCircuit(layers, is_adaptive, inModes, 
        outModes, ancModes, idlModes, num_of_modes)
        envModes = [m for m in 1:num_of_modes if m ∉ ancModes && m ∉ inModes]
        preps = SymplecticLayer(
            length(ancModes),
            x -> phaseShifting([m ∈ ancModes ? x[findfirst(==(m), ancModes)] : 0 for m in 1:num_of_modes]...)
        )
        postps = SymplecticLayer(
            length(idlModes),
            x -> phaseShifting([m ∈ idlModes ? x[findfirst(==(m), idlModes)] : 0 for m in 1:num_of_modes]...)
        )
        len_of_F = 2 * length(outModes) * length(idlModes)
        raw = *(layers...)
        A = SymplecticLayer(len_of_F, 
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
        num_of_args, f = final.num_of_args, final.f
        return new(layers, is_adaptive, inModes,
        outModes, ancModes, idlModes, envModes,
        num_of_args, num_of_modes, f)
    end
end

(sc::SymplecticCircuit)(x) = sc.f(x)

function symplecticLayer(ω::Float64, cr::CoupledResonators)
    return SymplecticLayer(
        0,
        x -> dilatedScatteringMatrix(ω, cr)
    )
end






