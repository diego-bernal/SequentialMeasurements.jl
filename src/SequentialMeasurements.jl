module SequentialMeasurements

using Reexport

# Load core types and modules in dependency order
include("MeasurementAnalysis.jl")  # Renamed from OutputProcessing
include("CoreTypes.jl")
include("BathOperators.jl")
include("ControlledEvolution.jl")
include("SpectralRepresentation.jl")
include("ParamDataBridge2.jl")
include("Initialization.jl")
include("QuantumProtocols.jl")
include("SimulationExecution.jl")

# Re-export all modules
@reexport using .MeasurementAnalysis
@reexport using .CoreTypes
@reexport using .BathOperators
@reexport using .ControlledEvolution
@reexport using .SpectralRepresentation
@reexport using .ParamDataBridge2
@reexport using .Initialization
@reexport using .QuantumProtocols
@reexport using .SimulationExecution

export initialize_noise, initialize_ctrl, initialize_measurement_scheme, run_simulation

end # module SequentialMeasurements
