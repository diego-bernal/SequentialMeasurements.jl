using QuantumOptics
using SequentialMeasurements

export error_virtualization_realization

"""
Realization function for error virtualization protocol.
This defines how a single realization of the protocol evolves.
"""
function error_virtualization_realization(
    params_operators::ParamsOperators, 
    control_params::CDDParams,
    ξs_t,
    ξb_t::Vector,
    measurement_scheme::MeasurementScheme
    )

    (; Ω, g, γs, γb, rhoq, rhob, sx, sy, sz, X, Y, Z) = params_operators
    (; n, t0, τg, nsteps, M_max) = control_params

    noise_operator_params = NoiseOperatorParams(Ω, g, γs, γb, X, Y, ξs_t, ξb_t)
    H = construct_hamiltonian(sz, noise_operator_params)

    # If no bath qubits, just evolve the system qubit
    if isempty(Ω)
        final_states = periodic_control(M_max, n, t0, τg, nsteps, rhoq, rhoq, H, [sx])
    else
        b = basis(X[1])
        final_states = periodic_control(M_max, n, t0, τg, nsteps, rhoq⊗rhob, rhoq, H, [sx⊗one(b)])
    end

    return process_measurement(measurement_scheme, final_states)
end