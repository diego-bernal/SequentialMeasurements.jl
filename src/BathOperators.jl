module BathOperators

using QuantumOpticsBase, QuantumOptics, StaticArrays, LinearAlgebra, Distributions

export bath_correlations, average_bath_correlations, NoiseOperatorParams, createNoiseOperator, sample_custom_pdf,
    construct_hamiltonian


struct NoiseOperatorParams
    Ω::Vector{Float64}
    g::Vector{Float64}
    γs::Float64              
    γb::Float64     
    X::Vector{Operator}
    Y::Vector{Operator}
    Z::Vector{Operator}
    ξs_t::Function              
    ξb_t::Vector{Function} 
end


function createNoiseOperator(params::NoiseOperatorParams)

    Nb = length(params.Ω)
    b = basis(params.X[1])

    return function(t)

        # Nb = 0 case
        if isempty(params.Ω)  
            return sqrt(params.γs)/2 * params.ξs_t(t)
        end

        B_t = sqrt(params.γs)/2 * params.ξs_t(t) * one(b)
        
        for i in 1:Nb
            ξbj_t = params.ξb_t[i](t)
            cos_term = cos(params.Ω[i]*t + sqrt(params.γb) * ξbj_t) * params.X[i]
            sin_term = -sin(params.Ω[i]*t + sqrt(params.γb) * ξbj_t) * params.Y[i]
            B_t += params.g[i] * (cos_term + sin_term)
        end
        
        return B_t
    end
end


# TODO: the current implementation is too specific to pure dephasing noise, 
# it should be more general. It also should be possible to give as an optional 
# argument to work on the interaction picture or not.
function construct_hamiltonian(sz::Operator, params::NoiseOperatorParams)

    # Nb = 0 case
    if isempty(params.Ω)  
        # Only system contribution with classical noise
        return sqrt(params.γs) / 2 * TimeDependentSum(params.ξs_t => sz)
    end
    
    b = basis(params.X[1])
    
    # system contribution
    H = sqrt(params.γs) / 2 * TimeDependentSum(params.ξs_t => sz ⊗ one(b))

    # loop to add terms for each bath qubit
    for j in 1:length(params.Ω)
        ξj = params.ξb_t[j]
        cos_func = t -> cos(params.Ω[j]*t + sqrt(params.γb)*ξj(t))
        sin_func = t -> -sin(params.Ω[j]*t + sqrt(params.γb)*ξj(t))

        H += params.g[j] * TimeDependentSum(cos_func => sz ⊗ params.X[j], sin_func => sz ⊗ params.Y[j])
    end

    return H
end


function bath_correlations(tspan, rhob, B)
    
    n = length(tspan)
    corr = zeros(ComplexF64, n, n)
    
    B_cache = [B(t) for t in tspan]
    
    for i in 1:n
        for j in i:n  # Only iterate over j >= i
            value = tr(B_cache[i] * B_cache[j] * rhob)
            corr[i, j] = value
            if i != j
                # use symmetry properties to fill the other half
                real_part = real(value)
                imag_part = -imag(value)
                corr[j, i] = complex(real_part, imag_part)
            end
        end
    end
    
    return tspan, tspan, corr
end


function average_bath_correlations(tspan, rhob, noise_operators)

    n = length(tspan)
    n_trajectories = length(noise_operators)
    total_corr = zeros(ComplexF64, n, n)

    for Bt in noise_operators
        _, _, corr = bath_correlations(tspan, rhob, Bt)
        total_corr .+= corr
    end

    avg_corr = total_corr ./ n_trajectories
    
    return tspan, tspan, avg_corr
end


function sample_custom_pdf(N, f; kwargs...)

    g = Normal(0, 0.3)
    c = 1000.0

    samples = []
    while length(samples) < N
        x = rand(g)
        u = rand(Uniform(0, c * pdf(g, x)))
        if u <= f(x; kwargs...)  # acceptance condition
            push!(samples, x)
        end
    end
    return abs.(samples)
end


end # module