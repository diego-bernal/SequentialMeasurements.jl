module CoreTypes

using QuantumOptics
using OrderedCollections

export ParamsOperators, NoiseParams, CDDParams, MeasurementScheme

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

"""
Measurement scheme parameters for processing simulation outputs
"""
struct MeasurementScheme
    fout::Union{Function, Nothing}  # Function to process states after each unitary block
    conditional_evolution::Bool      # Whether to use conditional evolution
    measurement_basis::Union{Operator, Nothing}  # Measurement basis for conditional evolution
    
    # Constructor with default values
    function MeasurementScheme(;
        fout::Union{Function, Nothing}=nothing,
        conditional_evolution::Bool=false,
        measurement_basis::Union{Operator, Nothing}=nothing
    )
        new(fout, conditional_evolution, measurement_basis)
    end
end

end # module CoreTypes