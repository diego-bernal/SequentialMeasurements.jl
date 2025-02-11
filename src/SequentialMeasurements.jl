module SequentialMeasurements

using Reexport

using QuantumOptics, StaticArrays, LinearAlgebra
using OrderedCollections
using Base.Threads, Logging
using Distributions # what functions am I using from this module?
import Random: randn
import Interpolations: cubic_spline_interpolation
import StatsBase: sample, ProbabilityWeights
import Statistics: mean


# Load custom modules
include("ParamDataBridge2.jl")
include("SpectralRepresentation.jl")
include("ControlledEvolution.jl")
include("BathOperators.jl")
@reexport using .ParamDataBridge2
@reexport using .SpectralRepresentation
@reexport using .ControlledEvolution
@reexport using .BathOperators


# Including lower-level functionalities
include("error_virtualization.jl")
export NoiseParams, initialize_simulation, multi_realization


# in the future I may want to consider putting the helper functions in a separate 
# module StochasticUtils.jl or StochasticCtrl.jl (StochasticControl.jl). A module 
# like this could contain SpectralRepresentation.jl, BathOperators.jl and 
# ControlledEvolution.
# Now I'm thinking of a separate module called SequentialMeasurements.jl

function initialize_ctrl(
    ctrl_type::Symbol, 
    ctrl_params_dict::Dict{Symbol, T} 
    ) where T

    if ctrl_type == :cdd
        return CDDParams(;ctrl_params_dict...)
    else
        error("Unknown simulation type: $ctrl_type")
    end

end


function initialize_noise(
    noise_type::Symbol, 
    noise_params_dict=Dict()
    ) 
    if noise_type == :lf_gaussian
        return NoiseParams(;noise_params_dict...)
    else
        error("Unknown simulation type: $noise_type")
    end
end


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


# High-level simulation functions
function run_simulation(
    num_realizations::Int64,  
    noise_params::NoiseParams,
    ctrl_params::CDDParams,
    data_collection::String;
    save_output=false
    )

    @info "\n\n--- New Simulation Initialization ---"
    @info "Preparing new simulation parameters and operators..."

    timed_result = @timed initialize_simulation(noise_params)
    params_operators = timed_result.value
    elapsed_time = timed_result.time

    @info "Initialization completed in $elapsed_time seconds."

    function serial_multi_realization(num_realizations, params_operators, noise_params, ctrl_params)
        start = time()
        average_fidelities = multi_realization(num_realizations, params_operators, noise_params, ctrl_params)
        mean_runtime = round((time() - start)/num_realizations, digits=3)
        @info "Thread $(threadid()) | Realizations run: $num_realizations | Mean realization runtime: $mean_runtime s"
        return average_fidelities
    end

    function distribute_tasks(num_realizations, n_cores)
        base = div(num_realizations, n_cores)
        remainder = num_realizations % n_cores
        return [base + (i <= remainder) for i in 1:n_cores]
    end

    n_cores = nthreads()
    if n_cores == 1
        @info "Running serial simulation for $num_realizations realizations."
        average_fidelities = serial_multi_realization(num_realizations, params_operators, noise_params, ctrl_params)
    else
        @info "Running parallel simulation with $n_cores/$(length(Sys.cpu_info())) cores for $num_realizations realizations."
        results = Vector{Vector{Float64}}()

        n_split = distribute_tasks(num_realizations, n_cores)
        @threads for n_runs in n_split
            local_result = serial_multi_realization(n_runs, params_operators, noise_params, ctrl_params)
            push!(results, local_result)
        end
        average_fidelities = mean(results)
    end

    if save_output
        params_dict = create_params_dict(num_realizations, params_operators, noise_params, ctrl_params)
        save_dataset(data_collection, params_dict, average_fidelities) 
    end

end


#------------------------------------------------------------------------------
# required function to run Ramsey spectroscopy simulations

# TODO: It would be good to have some output data requirements as input to the 
# run_simulation function rather than having separate functions for different
# types of simulations. This would make the code more modular and easier to
# maintain.


function run_simulation_ramsey(
    num_realizations::Int64,  
    noise_params::NoiseParams,
    ctrl_params::CDDParams,
    data_collection::String;
    save_output=false
    )

    @info "\n\n--- New Simulation Initialization ---"
    @info "Preparing new simulation parameters and operators..."

    timed_result = @timed initialize_simulation(noise_params)
    params_operators = timed_result.value
    elapsed_time = timed_result.time

    @info "Initialization completed in $elapsed_time seconds."

    function serial_multi_realization(num_realizations, params_operators, noise_params, ctrl_params)
        start = time()
        average_fidelities = multi_realization(num_realizations, params_operators, noise_params, ctrl_params)
        mean_runtime = round((time() - start)/num_realizations, digits=3)
        @info "Thread $(threadid()) | Realizations run: $num_realizations | Mean realization runtime: $mean_runtime s"
        return average_fidelities
    end

    function distribute_tasks(num_realizations, n_cores)
        base = div(num_realizations, n_cores)
        remainder = num_realizations % n_cores
        return [base + (i <= remainder) for i in 1:n_cores]
    end

    n_cores = nthreads()
    if n_cores == 1
        @info "Running serial simulation for $num_realizations realizations."
        average_fidelities = serial_multi_realization(num_realizations, params_operators, noise_params, ctrl_params)
    else
        @info "Running parallel simulation with $n_cores/$(length(Sys.cpu_info())) cores for $num_realizations realizations."
        results = Vector{Vector{Float64}}()

        n_split = distribute_tasks(num_realizations, n_cores)
        @threads for n_runs in n_split
            local_result = serial_multi_realization(n_runs, params_operators, noise_params, ctrl_params)
            push!(results, local_result)
        end
        average_fidelities = mean(results)
    end

    if save_output
        params_dict = create_params_dict(num_realizations, params_operators, noise_params, ctrl_params)
        save_dataset(data_collection, params_dict, average_fidelities) 
    end

end

#------------------------------------------------------------------------------


export initialize_noise, initialize_ctrl, run_simulation, run_simulation_ramsey, create_params_dict

end # module SequentialMeasurements
