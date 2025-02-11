struct ParamsOperators
    Ω::Vector{Float64}       # Rabi frequencies for the system
    g::Vector{Float64}       # Coupling strengths between system states and bath qubits
    γs::Float64              # system dephasing rate
    γb::Float64              # bath dephasing rate
    rhoq::Operator           # Density matrix for the quantum system
    rhob::Union{Operator,Int}# Density matrix for the bath (Int when Nb=0)
    sx::Operator             # Pauli X operator for the system
    sz::Operator             # Pauli Z operator for the system
    X::Vector{Operator}      # X operators corresponding to various noise sources
    Y::Vector{Operator}      # Y operators corresponding to various noise sources
end


function realization(
    params_operators::ParamsOperators, 
    control_params::CDDParams,
    ξs_t,
    ξb_t::Vector
    )

    (; Ω, g, γs, γb, rhoq, rhob, sx, sz, X, Y) = params_operators
    (; n, t0, τg, nsteps, M_max) = control_params

    # For Nb=0 case, ξb_t will be empty but that's fine since NoiseOperatorParams and construct_hamiltonian handle it
    noise_operator_params = NoiseOperatorParams(Ω, g, γs, γb, X, Y, ξs_t, ξb_t)
    H = construct_hamiltonian(sz, noise_operator_params)

    # If no bath qubits, just evolve the system qubit
    if isempty(Ω)
        final_states = periodic_control(M_max, n, t0, τg, nsteps, rhoq, rhoq, H, [sx])
    else
        b = basis(X[1])
        final_states = periodic_control(M_max, n, t0, τg, nsteps, rhoq⊗rhob, rhoq, H, [sx⊗one(b)])
    end

    fidelities = gate_fidelity.(final_states, Ref(rhoq))
    return fidelities
end


function multi_realization(
    n_runs::Int64, 
    params_operators::ParamsOperators, 
    noise_params::NoiseParams, 
    control_params::CDDParams
    )
    
    (; Nb, psd_function, psd_kwargs, ωu, Nh) = noise_params
    ξ_sys, ξ_bath = noise_generation_loop(n_runs, Nb, psd_function, ωu, Nh; psd_kwargs...)

    accumulated_fidelities = zeros(Float64, control_params.M_max+1)
    
    for i in 1:n_runs
        # For Nb=0, ξ_bath will be empty array, so use empty vector for that case
        bath_noise = isempty(ξ_bath) ? Function[] : ξ_bath[i]
        fidelities = realization(params_operators, control_params, ξ_sys[i], bath_noise)
        accumulated_fidelities .+= fidelities
    end

    average_fidelities = accumulated_fidelities / n_runs

    return average_fidelities
end


function init_operators(Nb)

    b_qubit = SpinBasis(1 // 2)

    # system operators
    sx = sigmax(b_qubit)
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
        rhob = 1
    elseif Nb == 1
        # Single bath qubit case
        X = [sigmax(b_qubit)]
        Y = [sigmay(b_qubit)]
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

    return rhoq, rhob, sx, sz, X, Y

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


function initialize_simulation(noise_params::NoiseParams)

    (; Nb, psd_function, psd_kwargs) = noise_params

    Ω, g, γs, γb = init_params(Nb, psd_function; psd_kwargs...)
    rhoq, rhob, sx, sz, X, Y = init_operators(Nb)

    return ParamsOperators(Ω, g, γs, γb, rhoq, rhob, sx, sz, X, Y)
end

