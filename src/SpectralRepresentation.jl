module SpectralRepresentation

using FFTW, SpecialFunctions, Interpolations, Random
using ..CoreTypes

export spectral_representation, interpolate_noise, noise_generation_loop, group_array,
    periodogram, fourier_transform, lorentzian, lf_gaussian, hf_gaussian, lorentzian_article

    
function lorentzian_article(omega; P0=1.0, omega_c=200.0)
    return (P0 / (π * omega_c)) / (1 + (omega/omega_c)^2)
end

function lorentzian(ω; σ=1.0, θ=1.0, ω₀=0.0)
    return σ^2 / (θ^2 + (ω-ω₀)^2) # FWHM = 2θ
end

function lf_gaussian(ω; s=1, Δ0=100, Γ0=0.1)
    return (1 / gamma(0.5 * (s + 1))) * (Δ0 / Γ0) * (abs(ω) / Γ0)^s * exp(-ω^2 / Γ0^2)
end

function hf_gaussian(ω; Δ1=0.1/20, Γ1=0.018, ω1=(0.99)*4π/2)
    return (1 / sqrt(π)) * (Δ1 / Γ1) * (exp(-((ω - ω1)^2) / Γ1^2) + exp(-((ω + ω1)^2) / Γ1^2))
end

integral_lf_gaussian(ω; s=1, Δ0=100, Γ0=0.1) = lf_gaussian(ω; s=s, Δ0=Δ0, Γ0=Γ0)./ω.^2;


"""
    spectral_representation(t_list, omega_m, num_realizations, psd_function; 
                          kwargs...)

Generate a set of realizations of a stochastic process given a power spectral 
density (PSD) function.

# Arguments
- `t_list::Vector{Float64}`: The time points at which the realizations are 
   sampled. This vector defines the temporal domain of the realizations.
- `omega_m::Vector{Float64}`: The angular frequencies at which the PSD function 
   is evaluated. It is assumed that these are equally spaced in frequency, are 
   all positive, and the first term is different than zero.
- `num_realizations::Int64`: The number of noise realizations to generate.
- `psd_function::Function`: A function that takes a frequency (and optionally 
   additional keyword arguments) and returns the PSD value at that frequency.
- `kwargs...`: Optional keyword arguments to pass to `psd_function`.

# Returns
- A matrix where each column is a realization of the stochastic process, sampled 
   at the time points specified in `t_list`.

# Notes
The function uses the spectral representation method to generate the realizations. 
This involves generating a set of random phases, constructing a complex signal in 
the frequency domain, and then taking the inverse Fourier transform to obtain the 
realizations in the time domain.

"""
function spectral_representation(
    t_list, 
    omega_m,
    num_realizations::Int64,  
    psd_function::Function;
    kwargs...)
    
    nh = length(omega_m)
    nt = length(t_list)
    Δω = omega_m[1]
    # Δt = t_list[2]-t_list[1]

    sqrt_Δω = sqrt(Δω / π)
    sqrt_2 = sqrt(2)
    exp_factor = 2π*im

    an_list = sqrt_Δω * sqrt.( (psd_function.(omega_m; kwargs...)))
    bn_matrix = zeros(ComplexF64, nt, num_realizations)

    random_phases = exp.(exp_factor * rand(Float64, (nh, num_realizations)))
    bn_matrix[2:(nh+1),:] .= sqrt_2 .* an_list .* random_phases

    fft_result = nt .* real.(ifft(bn_matrix, (1,)))

    return fft_result
end


"""
    interpolate_noise(time_vector::Vector{Float64}, noise_matrix::Matrix{Float64})

Generate interpolation functions for each noise realization.

# Arguments
- `time_vector`: A vector representing the time points.
- `noise_matrix`: A matrix where each column is a noise realization.

# Returns
- A vector of functions, each capable of interpolating the noise over time for a given realization.
"""
function interpolate_noise(time_vector, noise_matrix)
    num_realizations = size(noise_matrix, 2)
    inter_noise = Vector{Any}(undef, num_realizations)

    for i in 1:num_realizations
        inter_noise[i] = CubicSplineInterpolation(time_vector, noise_matrix[:, i], bc=Line(OnGrid()))
    end

    noise_functions = [t -> inter_noise[i](t) for i in 1:num_realizations]

    return noise_functions
end


"""
    group_array(original_array::Vector{T}, group_size::Int) -> Vector{Vector{T}}

Partition a one-dimensional array into a vector of vectors, each of size `group_size`.

# Arguments
- `original_array::Vector{T}`: The one-dimensional array to be partitioned, 
    where `T` is any type.
- `group_size::Int`: The size of each partition. The length of `original_array` 
    must be divisible by `group_size`.

# Returns
- A vector of vectors (`Vector{Vector{T}}`), where each inner vector is of length 
    `group_size`. The total number of inner vectors equals the length of 
    `original_array` divided by `group_size`.

# Example
```julia
original_array = [1, 2, 3, 4, 5, 6]
group_size = 2
grouped_arrays = group_array(original_array, group_size)
# Output: [[1, 2], [3, 4], [5, 6]]
````
"""
function group_array(original_array::Vector{T}, group_size::Int) where T
    # Handle empty array or group_size = 0 case
    if isempty(original_array) || group_size == 0
        return Vector{Vector{T}}()
    end

    # Calculate the number of groups based on the length of the original array and group_size
    num_groups = div(length(original_array), group_size)

    # Check divisibility to ensure a clean partition is possible
    if length(original_array) % group_size != 0
        error("The length of the original array is not divisible by group_size")
    end
    
    # Reshape the array into a 2D array of dimensions (group_size, num_groups)
    reshaped_array = reshape(original_array, group_size, num_groups)
    
    # Create a vector of vectors from the reshaped array
    grouped_arrays = [reshape(reshaped_array[:, i], group_size) for i in 1:num_groups]
    
    return grouped_arrays
end


"""
    noise_generation_loop(n_runs::Int, Nb::Int, psd_function::Function, 
                          ωu::Float64, Nh::Int; kwargs...)

Generate system and bath noise by performing a spectral representation of a 
power spectral density (PSD) function.

# Arguments
- `n_runs::Int`: The number of realizations of the noise process.
- `Nb::Int`: The number of bath qubits.
- `psd_function::Function`: The power spectral density function.
- `ωu::Float64`: The upper limit of the frequency sampling range.
- `Nh::Int`: The number of harmonics to consider.

# Keyword Arguments
- `kwargs`: Additional parameters to be passed to the `psd_function`.

# Returns
- `ξ_sys`: Array of system noise functions.
- `ξ_bath`: An array of arrays of length `Nb` of noise functions.

# Example
```julia
ξ_sys, ξ_bath = noise_generation_loop(1000, 2, psd_function, 2π*1.0, Int(2.0^12))
```
""" 
function noise_generation_loop(
    n_runs::Int, 
    Nb::Int,
    psd_function::Function,
    ωu::Float64, 
    Nh::Int;
    kwargs...)

    Δω = ωu/Nh
    T0 = 2π/Δω
    MM = 2*Nh
    ω = range(Δω, stop=ωu, length=Nh)
    t = range(0, stop=T0, length=MM)

    # System noise is always generated
    system_noise = spectral_representation(t, ω, n_runs, psd_function; kwargs...)
    ξ_sys = interpolate_noise(t, system_noise)

    # Bath noise only if Nb > 0
    if Nb == 0
        ξ_bath = Vector{Function}[]  # Return empty array of arrays for Nb=0
    else
        integral_psd_function(ω; kwargs...) = psd_function(ω; kwargs...)./ω.^2
        bath_noise = spectral_representation(t, ω, n_runs * Nb, integral_psd_function; kwargs...)
        ξ_bath = group_array(interpolate_noise(t, bath_noise), Nb)
    end

    return ξ_sys, ξ_bath
end


"""
    fourier_transform(t::Vector{Float64}, sampled_time_function::Vector{Float64})

Compute the Fourier transform of a sampled time function and return the angular 
frequencies and the normalized spectrum.

# Arguments
- `t::Vector{Float64}`: The time vector.
- `sampled_time_function::Vector{Float64}`: The values of the sampled time function.

# Returns
- `freqs::Vector{Float64}`: The angular frequencies.
- `spectrum::Vector{Float64}`: The normalized spectrum.

# Example
t = 0:Δt:tmax
sampled_function = sin.(2π * f * t)  # Example sampled function
freqs, spectrum = compute_fourier_transform(t, sampled_function)
"""
function fourier_transform(t::Vector{Float64}, sampled_time_function::Vector{Float64})
    # Ensure Δt is consistently calculated from t
    Δt = t[2] - t[1]  # Assuming uniform spacing
    # tmax = maximum(t)
    nt = length(t)

    # Sampling frequency
    fs = 1/Δt

    # Angular frequencies
    freqs = 2π * fftshift(fftfreq(nt, fs))

    # Compute the Fourier transform and normalize
    spectrum = Δt * abs.(fftshift(fft(sampled_time_function)))

    return freqs, spectrum
end



function periodogram(signal::Vector, fs::Float64)
    """
    arguments
    - `signal::Array`: An array containing the data
    - `fs::Float64`: Sampling rate (Hz). Number of points per unit time.

    returns
    - `w::Array{Float64}`: Fourier frequencies at which the periodogram is evaluated
    - `I_w::Array{Float64}`: The periodogram at frequences `w`
    """
    n = length(signal)
    w = 2π*rfftfreq(n, fs) # the 2π is necessary because I'm considering angular frequencies
    I_w = (1/fs) .* abs.(rfft(signal)).^2 ./ n # the periodogram has to be multiplied by Δt=1/fs

    return w, I_w
end


end #module