# Script that runs the error virtualization protocol simulation
using QuantumOptics
using SequentialMeasurements
using .SpectralRepresentation: lf_gaussian
using .MeasurementAnalysis: process_fidelity_with_plus

ctrl_params_dict = Dict(
    :n => 0,           # CDD order
    :t0 => 0.0,       # initial time
    :τg => 1.0,       # gate duration
    :nsteps => 10,    # number of time steps per free evolution region
    :M_max => 60      # number of control cycles
)

noise_params_dict = Dict(
    :Nb => 0,                             # number of bath qubits
    :psd_function => lf_gaussian,         # spectral density function
    :psd_kwargs => Dict(:s => 1.0, :Δ0 => 50, :Γ0 => 0.1),
    :ωu => 2π*1.0,                       # upper cutoff sampling frequency
    :Nh => Int(2.0^12)                   # number of frequncy sampling harmonics 
)

measurement_scheme_dict = Dict(
    :fout => process_fidelity_with_plus,  # From MeasurementAnalysis module
    :conditional_evolution => false
)

num_realizations = Int(1e1)
save_output = false
data_collection = @__DIR__

function main()
    control = initialize_ctrl(:cdd, ctrl_params_dict)
    measurement_scheme = initialize_measurement_scheme(measurement_scheme_dict)
    
    noise_params_sets = dict_list_extended(noise_params_dict)
    for noise_params in noise_params_sets
        noise = initialize_noise(:lf_gaussian, noise_params)
        run_simulation(
            num_realizations, 
            noise, 
            control, 
            data_collection,
            measurement_scheme,
            pure_dephasing_evolution;  # Using the generalized protocol directly
            save_output=save_output
        ) 
    end
end

main()

