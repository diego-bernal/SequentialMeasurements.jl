module SimulationExecution

using Base.Threads
using Logging
using Statistics: mean
using OrderedCollections
using ..CoreTypes
using ..MeasurementAnalysis
using ..Initialization
using ..ParamDataBridge2
using ..BathOperators
using ..SpectralRepresentation
using ..ControlledEvolution

export run_simulation, multi_realization, create_params_dict

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

function process_simulation_results(results, measurement_scheme::MeasurementScheme)
    if measurement_scheme.conditional_evolution
        # For conditional evolution, results is an array where each element is a tuple (states/fouts, outcomes)
        # Return raw results without averaging - one tuple per realization
        states_or_fouts = [r[1] for r in results]  # Array of arrays of states/fouts
        outcomes = [r[2] for r in results]         # Array of arrays of binary outcomes
        return states_or_fouts, outcomes
    else
        # For unconditional evolution, average across realizations
        if measurement_scheme.fout === nothing
            # Average states at each time point
            n_times = length(first(results))
            mean_states = [mean([results[r][t] for r in 1:length(results)]) for t in 1:n_times]
            return mean_states
        else
            # Average fout values at each time point
            n_times = length(first(results))
            mean_fouts = [mean([results[r][t] for r in 1:length(results)]) for t in 1:n_times]
            return mean_fouts
        end
    end
end

function multi_realization(
    n_runs::Int64, 
    params_operators::ParamsOperators, 
    noise_params::NoiseParams, 
    control_params::CDDParams,
    measurement_scheme::MeasurementScheme,
    realization_function::Function
    )
    
    (; Nb, psd_function, psd_kwargs, ωu, Nh) = noise_params
    ξ_sys, ξ_bath = noise_generation_loop(n_runs, Nb, psd_function, ωu, Nh; psd_kwargs...)
    
    results_list = Vector{Any}(undef, n_runs)
    
    for i in 1:n_runs
        bath_noise = isempty(ξ_bath) ? Function[] : ξ_bath[i]
        results_list[i] = realization_function(
            params_operators, 
            control_params,
            measurement_scheme,
            ξ_sys[i], 
            bath_noise
        )
    end

    return results_list
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
        partial_results = Vector{Vector}()

        n_split = distribute_tasks(num_realizations, n_cores)
        @threads for n_runs in n_split
            local_result = serial_multi_realization(n_runs, params_operators, noise_params, ctrl_params, measurement_scheme, realization_function)
            push!(partial_results, local_result)
        end
        results = vcat(partial_results...)
    end

    # Process results based on measurement scheme
    final_results = process_simulation_results(results, measurement_scheme)

    if save_output
        params_dict = create_params_dict(num_realizations, params_operators, noise_params, ctrl_params)
        save_dataset(data_collection, params_dict, final_results) 
    end

    return final_results
end

end # module SimulationExecution