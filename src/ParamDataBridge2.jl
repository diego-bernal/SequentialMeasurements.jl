module ParamDataBridge2


using Dates, JLD2, DrWatson, DataFrames, OrderedCollections


export dict_list_extended, save_dataset, load_dataset, sync_data!, delete_dataset!, delete_data_collection!


# extendend dict_list to handle nested dictionaries
function dict_list_extended(d::Dict)
    # Separate the entries into those that are dictionaries and those that are not
    dict_entries = Dict{Symbol,Any}()
    other_entries = Dict{Symbol,Any}()
    
    for (key, value) in d
        if isa(value, Dict)
            dict_entries[key] = value
        else
            other_entries[key] = value
        end
    end
    
    # Generate combinations for the non-dictionary entries
    other_combinations = dict_list(other_entries)
    
    # If there are no nested dictionaries, return the combinations from dict_list
    if isempty(dict_entries)
        return other_combinations
    end
    
    # Otherwise, recursively apply dict_list_modified to each nested dictionary
    nested_combinations = Dict{Symbol, Vector{Dict}}()
    for (key, value) in dict_entries
        nested_combinations[key] = dict_list_extended(value)
    end
    
    # Now, combine the results from other_combinations with nested_combinations
    final_combinations = Vector{Dict}()
    for combo in other_combinations
        # Start with a base combination
        temp_combinations = [deepcopy(combo)]
        for (key, values) in nested_combinations
            new_combos = []
            for val in values
                for existing in temp_combinations
                    push!(new_combos, merge(deepcopy(existing), Dict(key => val)))
                end
            end
            temp_combinations = new_combos
        end
        append!(final_combinations, temp_combinations)
    end
    
    return final_combinations
end


"""
    save_dataset(data_collection::String,
                 params_dict::AbstractDict{String, Any},
                 data::Any;
                 param_dir::Union{String, Nothing}=nothing, 
                 data_dir::String="dataset", 
                 base_filename::Union{String, Nothing}=nothing, 
                 time_stamp::Bool=true)

Save a dataset and associated parameters, organizing the files according to
the specified directories and options.

# Arguments
- `data_collection::String`: Path to the main directory where the data and
  parameters should be saved.
- `params_dict::AbstractDict{String, Any}`: Dictionary containing key-value
  pairs representing parameters associated with the dataset.
- `data::Any`: The actual data to be saved, which could be any Julia object.

# Keyword Arguments
- `param_dir::Union{String, Nothing}=nothing`: Directory name within
  `data_collection` where the parameters file will be saved. If `nothing`,
  the parameters file is saved directly within `data_collection`.
- `data_dir::String="dataset"`: Directory name within `data_collection`
  where the data file will be saved.
- `base_filename::Union{String, Nothing}=nothing`: Base name for the files.
  If `nothing`, files will only include the timestamp if `time_stamp` is
  `true`.
- `time_stamp::Bool=true`: Whether to include a timestamp in the filenames
  for uniqueness.

# Details
1. **Path Setup:** Combines the `data_collection` path with `param_dir` and
   `data_dir` to create the full paths where the files will be saved. Uses
   `mkpath` to ensure the required directories exist.
2. **File Naming:** Generates filenames using the optional `base_filename`
   and a timestamp (if `time_stamp` is `true`). Appends `_params.jld2` for
   the parameters file and `_data.jld2` for the data file.
3. **Saving the Data and Parameters:** Saves the `data` object using the
   `@save` macro. Updates `params_dict` to include the relative path to the
   data file. Uses `wsave` to write the parameters dictionary.

# Example
```julia
params = Dict("learning_rate" => 0.01, "batch_size" => 32)
data = rand(100, 10)

save_dataset("results/experiment1", params, data;
    param_dir="parameters",
    data_dir="data_files",
    base_filename="trial1",
    time_stamp=true)
"""
function save_dataset(data_collection::String,
    params_dict::AbstractDict{String, Any},
    data::Any; 
    param_dir::Union{String, Nothing}=nothing, 
    data_dir::String="dataset", 
    base_filename::Union{String, Nothing}=nothing, 
    time_stamp::Bool=true
    )

    param_full_path = isnothing(param_dir) ? data_collection : joinpath(data_collection, param_dir)
    data_full_path = joinpath(data_collection, data_dir)

    mkpath(param_full_path)
    mkpath(data_full_path)

    # timestamp for uniqueness, if enabled
    timestamp = time_stamp ? "_$(Dates.format(Dates.now(), "yyyymmdd_HHMMSS"))" : ""

    base_filename_part = isnothing(base_filename) ? "" : "$(base_filename)_"
    params_filename = "$(base_filename_part)$(timestamp)_params.jld2"
    data_filename = "$(base_filename_part)$(timestamp)_data.jld2"

    params_path_filename = joinpath(param_full_path, params_filename)
    data_path_filename = joinpath(data_full_path, data_filename)

    @save data_path_filename data

    params_dict["data_path"] = joinpath(data_dir, data_filename)
    wsave(params_path_filename, params_dict)

end


"""
    load_dataset(df::DataFrame, row::Int64)

Load a dataset by combining the directory path with a relative data path
found in the specified DataFrame row.

# Arguments
- `df::DataFrame`: The DataFrame containing path information in specific
  columns.
- `row::Int64`: The row number to read from the DataFrame.

# Details
1. **Row Validation:** Verifies that the specified `row` is within the valid
   range of the DataFrame.
2. **Extract Paths:** Retrieves the `:path` column (full directory path) and
   the `:data_path` column (relative data path) from the specified row.
3. **Validate Paths:** Ensures both paths are valid strings and are not
   missing or null.
4. **Combine Paths:** Extracts the directory part of the `:path` value using
   `dirname` and combines it with the relative `:data_path`.
5. **Load Dataset:** Loads the dataset from the combined path using
   DrWatson's `wload` function.

# Example
```julia
# Example DataFrame setup
using DataFrames
df = DataFrame(path=["/full/directory/path/file1.jld2",
                     "/full/directory/path/file2.jld2"],
               data_path=["relative/file1.jld2", "relative/file2.jld2"])

# Load a specific dataset using the load_dataset function
dataset = load_dataset(df, 1)
"""
function load_dataset(df::DataFrame, row::Int64)

    if row < 1 || row > nrow(df)
        error("The specified row is out of bounds.")
    end

    full_path = df[row, :path]
    data_path = df[row, :data_path] 

    if typeof(full_path) != String || ismissing(full_path)
        error("Invalid or missing directory path in the specified row.")
    end

    if typeof(data_path) != String || ismissing(data_path)
        error("Invalid or missing data path in the specified row.")
    end

    directory = dirname(full_path)
    combined_path = joinpath(directory, data_path)

    # Load the dataset using DrWatson's wload function
    dataset = wload(combined_path)

    return dataset
end


# TODO I should think of a way to include the black list here by reading the DataFrame df.
function sync_data!(df::DataFrame)

    # Assume all 'path' entries share the same root directory for 'data_collection'
    if !isempty(df.path)
        first_path = df[1, "path"]
        data_collection = dirname(first_path)  # Infer the collection folder from the first entry
        
        # Perform synchronization of data
        for row in eachrow(df)
            param_exists = isfile(row["path"])
            data_exists = isfile(row["data_path"])
            
            if param_exists && !data_exists
                rm(row["path"], force = true)
            elseif !param_exists && data_exists
                rm(row["data_path"], force = true)
            end
        end

        # Collect all files under data_collection and its subdirectories
        all_files = list_files_recursively(data_collection)

        # Get a list of all referenced files in the DataFrame
        referenced_files = unique(vcat(df.path, df.data_path))

        # Filter out the files that are not referenced and delete them
        for file in all_files
            if !(file in referenced_files)
                rm(file, force = true)
            end
        end

        results_file = joinpath(dirname(data_collection), "results_$(basename(data_collection)).jld2")
        if isfile(results_file)
            rm(results_file, force = true)
            collect_results!(data_collection)
        else
            collect_results(data_collection)
        end

    else
        println("DataFrame is empty or 'path' column is missing entries.")
    end

end


function list_files_recursively(directory::String)

    files = String[]
    for entry in readdir(directory, join=true)  # join=true to get full paths
        if isdir(entry)
            append!(files, list_files_recursively(entry))  # Recursive call for subdirectories
        else
            push!(files, entry)
        end
    end
    return files
end


function delete_dataset!(df::DataFrame, row::Int)

    # Validate the row index
    if row < 1 || row > nrow(df)
        error("Row index out of bounds.")
    end

    # Extract file paths from the DataFrame
    path_to_delete = df[row, "path"]
    
    # Delete the file in the 'path' column
    if isfile(path_to_delete)
        rm(path_to_delete, force=true)
        sync_data!(df)
    else
        sync_data!(df)
    end

end


function delete_data_collection!(data_collection::String)
    # Check if the directory exists
    if isdir(data_collection)
        # Recursively delete all contents of the directory
        for entry in readdir(data_collection, join=true)  # join=true to get full paths
            if isdir(entry)
                delete_data_collection!(entry)  # Recursive call for subdirectories
            else
                rm(entry, force=true)  # Remove files with force to avoid errors if the file is already in use
            end
        end
        # After clearing contents, remove the directory itself
        rm(data_collection, force=true)
        println("Deleted directory: ", data_collection)
    else
        println("Data collection does not exist: ", data_collection)
    end
end


end #module