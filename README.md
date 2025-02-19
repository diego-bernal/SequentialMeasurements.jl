# SequentialMeasurements.jl

[![Build Status](https://github.com/diego-bernal/SequentialMeasurements.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/diego-bernal/SequentialMeasurements.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/diego-bernal/SequentialMeasurements.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/diego-bernal/SequentialMeasurements.jl)

A Julia package for simulating quantum systems under sequential measurements with dynamical decoupling control, particularly focusing on quantum error virtualization protocols.

## Features

- Simulation of quantum systems with configurable number of bath qubits
- Support for different noise models through customizable power spectral density functions
- Implementation of Concatenated Dynamical Decoupling (CDD) sequences
- Flexible measurement schemes (fidelity measurements, expectation values)
- Parallel execution support for multiple realizations
- Data management and persistence capabilities

## Installation

To install the package for development:

```julia
julia> using Pkg
julia> Pkg.develop(url="https://github.com/diego-bernal/SequentialMeasurements.jl.git")
```

## Usage Example

Here's a basic example of running an error virtualization simulation:

```julia
using SequentialMeasurements

# Define control parameters
ctrl_params_dict = Dict(
    :n => 0,           # CDD order
    :t0 => 0.0,       # initial time
    :τg => 1.0,       # gate duration
    :nsteps => 10,    # time steps per evolution
    :M_max => 60      # control cycles
)

# Define noise parameters
noise_params_dict = Dict(
    :Nb => 0,                             # number of bath qubits (0 for classical noise)
    :psd_function => lf_gaussian,         # spectral density function
    :psd_kwargs => Dict(:s => 1.0, :Δ0 => 50, :Γ0 => 0.1),
    :ωu => 2π*1.0,                       # upper frequency cutoff
    :Nh => Int(2.0^12)                   # number of harmonics
)

# Define measurement scheme
measurement_scheme_dict = Dict(
    :fout => process_fidelity_with_plus,  # measurement processing function
    :conditional_evolution => false        # whether to use conditional evolution
)

# Initialize components
control = initialize_ctrl(:cdd, ctrl_params_dict)
noise = initialize_noise(:lf_gaussian, noise_params_dict)
measurement_scheme = initialize_measurement_scheme(measurement_scheme_dict)

# Run simulation
results = run_simulation(
    1000,                                 # number of realizations
    noise,
    control,
    "output_directory",                   # data collection path
    measurement_scheme,
    error_virtualization_realization;     # realization function 
    save_output=false                     # whether to save results to disk
)
```

### Available Measurement Options

The package provides several measurement processing functions:

- `process_fidelity_with_plus`: Measures fidelity with respect to the |+⟩ state
- `process_fidelity`: Generic fidelity measurement with respect to any target state
- `process_expectation`: Measures expectation value of an operator

Measurements can be performed in two modes:
- Unconditional evolution (`conditional_evolution = false`): Average over all possible measurement outcomes
- Conditional evolution (`conditional_evolution = true`): Perform measurements and update states based on outcomes

### Noise Models

The package includes several built-in spectral density functions:
- `lf_gaussian`: Low-frequency Gaussian noise
- `hf_gaussian`: High-frequency Gaussian noise
- `lorentzian`: Lorentzian spectral density
- `lorentzian_article`: Alternative Lorentzian parameterization

You can also implement custom spectral density functions for specific noise models.

## Contributing

Contributions are welcome! Here's how you can help:

1. Fork the repository
2. Create a new branch (`git checkout -b feature/your-feature`)
3. Make your changes
4. Run the tests (`julia> ]test SequentialMeasurements`)
5. Commit your changes (`git commit -am 'Add some feature'`)
6. Push to the branch (`git push origin feature/your-feature`)
7. Create a new Pull Request

### Development Setup

To set up the package for development:

1. Clone the repository:
```bash
git clone https://github.com/diego-bernal/SequentialMeasurements.jl.git
cd SequentialMeasurements.jl
```

2. Start Julia and enter package mode by pressing `]`, then:
```julia
pkg> activate .
pkg> instantiate
```

3. Run the tests:
```julia
pkg> test
```

### Adding New Features

- **Measurement Schemes**: Add new measurement schemes by creating a subtype of `MeasurementScheme` and implementing `process_measurement`
- **Noise Models**: Add new noise models by implementing new spectral density functions
- **Control Sequences**: Extend the control capabilities by adding new control sequence generators

## License

This project is licensed under the MIT License - see the LICENSE file for details.
