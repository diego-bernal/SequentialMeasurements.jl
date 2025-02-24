using QuantumOptics
using SequentialMeasurements

export error_virtualization_realization

"""
Realization function for error virtualization protocol.
This defines how a single realization of the protocol evolves, returning the final states.
"""
function error_virtualization_realization(
    params_operators::ParamsOperators, 
    control_params::CDDParams,
    measurement_scheme::MeasurementScheme,
    ξs_t,
    ξb_t::Vector
    )

    (; Ω, g, γs, γb, rhoq, rhob, sx, sy, sz, X, Y, Z) = params_operators
    (; n, t0, τg, nsteps, M_max) = control_params

    noise_operator_params = NoiseOperatorParams(Ω, g, γs, γb, X, Y, Z, ξs_t, ξb_t)
    H = construct_hamiltonian(sz, noise_operator_params)

    # If no bath qubits, just evolve the system qubit
    if isempty(Ω)
        return periodic_control(M_max, n, t0, τg, nsteps, rhoq, rhoq, H, [sx], measurement_scheme)
    else
        b = basis(X[1])
        return periodic_control(M_max, n, t0, τg, nsteps, rhoq⊗rhob, rhoq, H, [sx⊗one(b)], measurement_scheme)
    end
    
end