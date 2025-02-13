module MeasurementSchemes

using ..CoreTypes
using QuantumOptics

export accumulate_measurements, initialize_basic_measurements, fidelity_plus

"""
Generic function to accumulate and average measurements over multiple realizations.
"""
function accumulate_measurements(results_list)
    accumulated = zeros(Float64, length(first(results_list)))
    for results in results_list
        accumulated .+= results
    end
    return accumulated / length(results_list)
end

"""
Initialize basic measurement schemes that can be used across different simulations.
Returns common measurement operators (σx, σy, σz) as ExpectationMeasurement schemes.
"""
function initialize_basic_measurements()
    return init_basic_measurements()  # Call the implementation from CoreTypes
end

"""
Create a measurement scheme for fidelity with respect to the plus state.
This is commonly used in error virtualization protocols.
"""
function fidelity_plus()
    b_qubit = SpinBasis(1 // 2)
    plus = spinup(b_qubit) + spindown(b_qubit)
    normalize!(plus)
    return FidelityMeasurement(dm(plus))
end

end # module MeasurementSchemes