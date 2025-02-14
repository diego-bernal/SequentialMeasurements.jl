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
ctrl_params = Dict(
    :n => 0,           # CDD order
    :t0 => 0.0,       # initial time
    :τg => 1.0,       # gate duration
    :nsteps => 10,    # time steps per evolution
    :M_max => 60      # control cycles
)

# Define noise parameters
noise_params = Dict(
    :Nb => 2,                             # number of bath qubits
    :psd_function => lf_gaussian,         # spectral density function
    :psd_kwargs => Dict(:s => 1.0, :Δ0 => 50, :Γ0 => 0.1),
    :ωu => 2π*1.0,                       # upper frequency cutoff
    :Nh => Int(2.0^12)                   # number of harmonics
)

# Initialize control and measurement scheme
control = initialize_ctrl(:cdd, ctrl_params)
measurement_scheme::MeasurementScheme = FidelityWithPlus()  # Type-safe measurement scheme

# Run simulation
results = run_simulation(
    1000,                  # number of realizations
    initialize_noise(:lf_gaussian, noise_params),
    control,
    "output_directory",    # data collection path
    measurement_scheme,
    realization_function   # realization function specific to your protocol
)
```

### Available Measurement Schemes

The package provides several measurement schemes:

- `FidelityWithPlus()`: Measures fidelity with respect to the plus state (|+⟩)
- `FidelityMeasurement(target_state)`: Measures fidelity with respect to any target state
- `ExpectationMeasurement(operator)`: Measures expectation value of an operator

You can also use `initialize_basic_measurements()` to get common measurement operators (σx, σy, σz).

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
