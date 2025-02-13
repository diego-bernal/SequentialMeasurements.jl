module Initialization

using ..CoreTypes
using QuantumOptics
using Distributions
using ..BathOperators

export initialize_ctrl, initialize_noise, initialize_simulation, init_params, init_operators

function init_operators(Nb)
    b_qubit = SpinBasis(1 // 2)

    # system operators
    sx = sigmax(b_qubit)
    sy = sigmay(b_qubit)
    sz = sigmaz(b_qubit)

    plus = spinup(b_qubit) + spindown(b_qubit)
    normalize!(plus)

    rhoq = dm(plus) # system qubit initial condition
    normalize!(rhoq)

    # bath operators
    if Nb == 0
        # No bath case - only classical noise
        X = Operator[]
        Y = Operator[]
        Z = Operator[]
        rhob = 1
    elseif Nb == 1
        # Single bath qubit case
        X = [sigmax(b_qubit)]
        Y = [sigmay(b_qubit)]
        Z = [sigmaz(b_qubit)]
        rhob = thermalstate(sigmaz(b_qubit), 1e0)
    else
        # Multiple bath qubits case
        b_bath = tensor([b_qubit for _ = 1:Nb]...)
        x(i) = embed(b_bath, i, sigmax(b_qubit))
        y(i) = embed(b_bath, i, sigmay(b_qubit))
        z(i) = embed(b_bath, i, sigmaz(b_qubit))

        X = x.(1:Nb)
        Y = y.(1:Nb)
        Z = z.(1:Nb)
        
        rhob = tensor([thermalstate(sigmaz(b_qubit), 1e0) for _ = 1:Nb]...)
    end

    if Nb > 0
        normalize!(rhob)
    end

    return rhoq, rhob, sx, sy, sz, X, Y, Z
end

function init_params(Nb, psd_function; kwargs...)
    if Nb == 0
        Ω = Float64[]
        g = Float64[]
    else
        Ω = sample_custom_pdf(Nb, psd_function; kwargs...)
        g = rand(Uniform(0, 1), Nb)
    end
    
    γs = 0.01
    γb = 0.1  

    return Ω, g, γs, γb
end

function initialize_ctrl(
    ctrl_type::Symbol, 
    ctrl_params_dict::Dict{Symbol, T} 
    ) where T

    if ctrl_type == :cdd
        return CDDParams(;ctrl_params_dict...)
    else
        error("Unknown simulation type: $ctrl_type")
    end
end

function initialize_noise(
    noise_type::Symbol, 
    noise_params_dict=Dict()
    ) 
    if noise_type == :lf_gaussian
        return NoiseParams(;noise_params_dict...)
    else
        error("Unknown simulation type: $noise_type")
    end
end

function initialize_simulation(noise_params::NoiseParams)
    (; Nb, psd_function, psd_kwargs) = noise_params

    Ω, g, γs, γb = init_params(Nb, psd_function; psd_kwargs...)
    rhoq, rhob, sx, sy, sz, X, Y, Z = init_operators(Nb)

    return ParamsOperators(Ω, g, γs, γb, rhoq, rhob, sx, sy, sz, X, Y, Z)
end

end # module Initialization