# Script that receives the parameters and run the `main()` function which processes and save the results
using DrWatson
@quickactivate :ErrorVirtualization

ctrl_params_dict = Dict(
    :n => 0, #1,        # CDD order
    :t0 => 0.0,     # initial time
    :τg => 1.0,     # gate duration
    :nsteps => 10,  # number of time steps per free evolution region
    :M_max => 60 #60    # number of control cycles
)

noise_params_dict = Dict(
    :Nb => 0, #2, #6,                        # number of bath qubits
    :psd_function => lf_gaussian,        # spectral density function
    :psd_kwargs => Dict(:s => 1.0, :Δ0 => 50, :Γ0 => 0.1),
    # :psd_kwargs => Dict(:s => 1.0, :Δ0 => [50, 200], :Γ0 => [0.1, 0.2]),
    :ωu => 2π*1.0,                       # upper cutoff sampling frequency
    :Nh => Int(2.0^12)                   # number of frequncy sampling harmonics 
)

save_output = false #true
data_collection = datadir("sims", "local", "periodic_control", "lf_noise")

num_realizations = Int(1e1) # 2e3

function main()
    
    control = initialize_ctrl(:cdd, ctrl_params_dict)

    noise_params_sets = dict_list_extended(noise_params_dict)
    for noise_params in noise_params_sets
        noise = initialize_noise(:lf_gaussian, noise_params)
        run_simulation(num_realizations, noise, control, data_collection; save_output=save_output) 
    end

end

main()

