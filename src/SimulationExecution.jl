module SimulationExecution

using Base.Threads
using Logging
using Statistics: mean
using ..CoreTypes
using ..Initialization
using ..ParamDataBridge2
using ..BathOperators
using ..SpectralRepresentation
using ..MeasurementSchemes

export run_simulation, multi_realization

function multi_realization(
    n_runs::Int64, 
    params_operators::ParamsOperators, 
    noise_params::NoiseParams, 
    control_params::CDDParams,
    measurement_scheme::MeasurementScheme,
    realization_function::Function  # Generic realization function passed as parameter
    )
    
    (; Nb, psd_function, psd_kwargs, ωu, Nh) = noise_params
    ξ_sys, ξ_bath = noise_generation_loop(n_runs, Nb, psd_function, ωu, Nh; psd_kwargs...)
    
    results_list = Vector{Vector{Float64}}(undef, n_runs)
    
    for i in 1:n_runs
        bath_noise = isempty(ξ_bath) ? Function[] : ξ_bath[i]
        results_list[i] = realization_function(params_operators, control_params, ξ_sys[i], bath_noise, measurement_scheme)
    end

    return accumulate_measurements(results_list)
end

function serial_multi_realization(num_realizations, params_operators, noise_params, ctrl_params, measurement_scheme, realization_function)
    start = time()
    results = multi_realization(num_realizations, params_operators, noise_params, ctrl_params, measurement_scheme, realization_function)
    mean_runtime = round((time() - start)/num_realizations, digits=3)
    @info "Thread $(threadid()) | Realizations run: $num_realizations | Mean realization runtime: $mean_runtime s"
    return results
end

function distribute_tasks(num_realizations, n_cores)
    base = div(num_realizations, n_cores)
    remainder = num_realizations % n_cores
    return [base + (i <= remainder) for i in 1:n_cores]
end

function run_simulation(
    num_realizations::Int64,  
    noise_params::NoiseParams,
    ctrl_params::CDDParams,
    data_collection::String,
    measurement_scheme::MeasurementScheme,
    realization_function::Function;
    save_output=false
    )

    @info "\n\n--- New Simulation Initialization ---"
    @info "Preparing new simulation parameters and operators..."

    timed_result = @timed initialize_simulation(noise_params)
    params_operators = timed_result.value
    elapsed_time = timed_result.time

    @info "Initialization completed in $elapsed_time seconds."

    n_cores = nthreads()
    if n_cores == 1
        @info "Running serial simulation for $num_realizations realizations."
        results = serial_multi_realization(num_realizations, params_operators, noise_params, ctrl_params, measurement_scheme, realization_function)
    else
        @info "Running parallel simulation with $n_cores/$(length(Sys.cpu_info())) cores for $num_realizations realizations."
        partial_results = Vector{Vector{Float64}}()

        n_split = distribute_tasks(num_realizations, n_cores)
        @threads for n_runs in n_split
            local_result = serial_multi_realization(n_runs, params_operators, noise_params, ctrl_params, measurement_scheme, realization_function)
            push!(partial_results, local_result)
        end
        results = mean(partial_results)
    end

    if save_output
        params_dict = create_params_dict(num_realizations, params_operators, noise_params, ctrl_params)
        save_dataset(data_collection, params_dict, results) 
    end

    return results
end

end # module SimulationExecution