module QuantumProtocols

using QuantumOptics
using ..CoreTypes
using ..BathOperators
using ..ControlledEvolution

export pure_dephasing_evolution

"""
    pure_dephasing_evolution(params_operators, control_params, measurement_scheme, ξs_t, ξb_t)

General realization function for pure dephasing evolution protocols.
This defines how a single realization of a pure dephasing protocol evolves, returning the final states.

# Arguments
- `params_operators::ParamsOperators`: Contains system and bath operators and parameters
- `control_params::CDDParams`: Control sequence parameters
- `measurement_scheme::MeasurementScheme`: Configuration for measurements
- `ξs_t::Function`: System noise function
- `ξb_t::Vector{Function}`: Bath noise functions

# Returns
- Evolution results according to the measurement scheme configuration
"""
function pure_dephasing_evolution(
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

end # module