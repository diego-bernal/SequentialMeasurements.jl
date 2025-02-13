module CoreTypes

using QuantumOptics
using OrderedCollections

export ParamsOperators, MeasurementScheme, FidelityMeasurement, ExpectationMeasurement, 
       NoiseParams, CDDParams, create_params_dict, init_basic_measurements, process_measurement

struct ParamsOperators
    Ω::Vector{Float64}       # Rabi frequencies for the system
    g::Vector{Float64}       # Coupling strengths between system states and bath qubits
    γs::Float64             # system dephasing rate
    γb::Float64             # bath dephasing rate
    rhoq::Operator          # Density matrix for the quantum system
    rhob::Union{Operator,Int}# Density matrix for the bath (Int when Nb=0)
    sx::Operator            # Pauli X operator for the system
    sy::Operator            # Pauli Y operator for the system
    sz::Operator            # Pauli Z operator for the system
    X::Vector{Operator}     # X operators corresponding to various noise sources
    Y::Vector{Operator}     # Y operators corresponding to various noise sources
    Z::Vector{Operator}     # Z operators corresponding to various noise sources
end

# Abstract type for different measurement schemes
abstract type MeasurementScheme end

# Concrete implementation for fidelity measurements
struct FidelityMeasurement <: MeasurementScheme
    target_state::Operator
end

# Concrete implementation for expectation value measurements
struct ExpectationMeasurement <: MeasurementScheme
    operator::Operator
end

# Function to process measurements based on scheme
function process_measurement(scheme::FidelityMeasurement, final_states)
    return gate_fidelity.(final_states, Ref(scheme.target_state))
end

function process_measurement(scheme::ExpectationMeasurement, final_states)
    return real.(expect.(Ref(scheme.operator), final_states))
end

"""
Initialize basic measurement schemes that are commonly used across different simulations
"""
function init_basic_measurements()
    b_qubit = SpinBasis(1 // 2)
    
    # Common measurement operators
    sz_scheme = ExpectationMeasurement(sigmaz(b_qubit))
    sx_scheme = ExpectationMeasurement(sigmax(b_qubit))
    sy_scheme = ExpectationMeasurement(sigmay(b_qubit))
    
    return (sz=sz_scheme, sx=sx_scheme, sy=sy_scheme)
end

@kwdef struct NoiseParams
    Nb::Int64                         # number of bath qubits
    psd_function::Function            # spectral density function
    psd_kwargs::Dict{Symbol, Any}
    ωu::Float64                       # upper cutoff sampling frequency
    Nh::Int64                        # number of harmonics to sample
end

@kwdef struct CDDParams
    n::Int64                         # CDD order
    t0::Float64                      # initial time
    τg::Float64                      # gate duration
    nsteps::Int64                   # number of time steps per free evolution region
    M_max::Int64                    # number of control cycles
end

# Function to create a parameter dictionary for data storage
function create_params_dict(
    num_realizations::Int64,
    params_operators::ParamsOperators,
    noise_params::NoiseParams,
    control_params::CDDParams
    )

    result_dict = OrderedDict{String, Any}()
    result_dict["num_realizations"] = num_realizations

    function add_selected_fields_to_dict(data_struct, dict, selected_fields)
        for field in selected_fields
            if hasfield(typeof(data_struct), field)
                dict[string(field)] = getfield(data_struct, field)
            end
        end
    end

    function add_fields_to_dict(data_struct, dict)
        for field in fieldnames(typeof(data_struct))
            value = getfield(data_struct, field)
            key = string(field)
            if value isa Function
                dict[key] = string(value)
            elseif value isa AbstractDict
                # Flatten the dictionary: append sub-key to main key
                for (k, v) in value
                    dict["$key.$k"] = v
                end
            else
                dict[key] = value
            end
        end
    end

    selected_operator_fields = [:Ω, :g, :γs, :γb]
    add_selected_fields_to_dict(params_operators, result_dict, selected_operator_fields)
    add_fields_to_dict(control_params, result_dict)
    add_fields_to_dict(noise_params, result_dict)

    return result_dict
end

end # module CoreTypes