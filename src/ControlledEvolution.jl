module ControlledEvolution

using QuantumOpticsBase, QuantumOptics, StaticArrays, LinearAlgebra, Distributions, Random
using ..CoreTypes

export cdd, cdd_filter_function, reset_qubit, cdd_final_state, cdd_time_evolution, 
       periodic_control_final_state, periodic_control_time_evolution, periodic_control,
       periodic_control_conditional, periodic_control_unconditional


"""
        cdd(n::Int)

Generates the `n`-th order concatenated dynamical decoupling (CDD) sequence.

# Arguments
- `n::Int`: Order of the CDD sequence.

# Returns
- An array representing the CDD sequence. Each element in the array represents a 
    pulse or a free evolution period. A `0` represents a free evolution period, and 
    a non-zero number represents  an instantaneous pulse raised to that power.

# Example
```julia
println(cdd(0))  # Prints: [0]
println(cdd(1))  # Prints: [0, 1, 0, 1]
println(cdd(2))  # Prints: [0, 1, 0, 2, 0, 1, 0, 2]
println(count(x -> x == 0, cdd(10))) # Prints the number of free evolution regions (2^n)
```
"""
function cdd(n::Int)
    
    if n == 0
        return [0]
    else
        prev = cdd(n - 1)  # Recursive call to get the previous array
        
        expanded = Vector{Int}(undef, 0)
        for val in prev
            if val == 0
                append!(expanded, [0, 1, 0, 1])  # Replace 0 with [0, 1, 0, 1]
            else
                push!(expanded, val) # Keep non-zero values as is
            end
        end

        # Process the expanded array to add adjacent non-zero numbers
        result = Vector{Int}()
        i = 1
        while i <= length(expanded)
            if i < length(expanded) && expanded[i] != 0 && expanded[i + 1] != 0
                # Add adjacent non-zero numbers
                push!(result, expanded[i] + expanded[i + 1])
                i += 2  # Skip the next number since it has been added
            else
                push!(result, expanded[i])
                i += 1
            end
        end

        return result
    end
end


function cdd_filter_function(ω, t, n)

    squaredTerm = 2^(2*n + 1) *(2/ω^2)
    sineSquaredTerm = sin(ω*t / 2^(n + 1))^2
    
    productSeries = 1.0 
    for j = 1:n
        productSeries *= sin(ω*t / 2^(j + 1))^2
    end

    result = squaredTerm * sineSquaredTerm * productSeries
    return result
end


function reset_qubit(ρsb::Operator, qubit_state::Operator)
    # If dimensions match, there's no bath to trace out
    if size(ρsb) == size(qubit_state)
        return qubit_state
    end
    return tensor(qubit_state, ptrace(ρsb, 1))
end


function cdd_final_state(
    n::Int64,
    t0::Float64,
    τg::Float64,
    nsteps::Int64, # number of steps for the free-evolution regions. Total number of steps = 2^n * nsteps
    rho0::Operator,
    H::TimeDependentSum,
    G::Union{SVector,Vector} # CDD basic operators
    )

    seq = cdd(n)
    τ0 = τg/2^n
    Δt = τ0/nsteps

    # t_collected = []
    # rho_collected = [] 

    idx = 1
    for val in seq
        if val == 0
            _, rho = timeevolution.master_dynamic([t0+(idx-1)*τ0:Δt:t0+idx*τ0-Δt;], rho0, H, [])
            # append!(t_collected, tout)
            # append!(rho_collected, rho)
            rho0 = last(rho)
            normalize!(rho0)
            idx += 1
        else
            rho0 = (G[1]^val)*rho0*dagger(G[1]^val) # implementation of dynamical decoupling gate
        end
    end

    return rho0
end


function cdd_time_evolution(
    n::Int64,
    t0::Float64,
    τg::Float64,
    nsteps::Int64, # number of steps for the free-evolution regions. Total number of steps = 2^n * nsteps
    rho0::Operator,
    H::TimeDependentSum,
    G::Union{SVector,Vector} # CDD basic operators
    )

    seq = cdd(n)
    τ0 = τg/2^n
    Δt = τ0/nsteps

    t_collected = []
    rho_collected = [] 

    idx = 1
    for val in seq
        if val == 0
            tout, rho = timeevolution.master_dynamic([t0+(idx-1)*τ0:Δt:t0+idx*τ0-Δt;], rho0, H, [])
            append!(t_collected, tout)
            append!(rho_collected, rho)

            rho0 = last(rho)
            normalize!(rho0)
            idx += 1
        else
            rho0 = (G[1]^val)*rho0*dagger(G[1]^val) # implementation of dynamical decoupling gate
        end
    end

    return t_collected, rho_collected
end


function periodic_control_final_state(
    M::Int64,
    n::Int64,
    t0::Float64,
    τg::Float64,
    nsteps::Int64, 
    rho0::Operator,
    rho_reset::Operator,
    H::TimeDependentSum,
    G::Union{SVector,Vector}
    )

    for m in 1:M
        rho0 = cdd_final_state(n, t0+(m-1)*τg, τg, nsteps, rho0, H, G)
    end

    rho0 = reset_qubit(rho0, rho_reset)
    
    rho0 = cdd_final_state(n, t0+M*τg, τg, nsteps, rho0, H, G)

    return rho0
end


function periodic_control_time_evolution(
    M::Int64,
    n::Int64,
    t0::Float64,
    τg::Float64,
    nsteps::Int64, 
    rho0::Operator,
    rho_reset::Operator,
    H::TimeDependentSum,
    G::Union{SVector,Vector}
)
    t_collected = Float64[]  # Use Float64 for consistency in time values
    rho_collected = Operator[]  # Assuming Operator is the correct type for rho

    for m in 1:M
        t, rho = cdd_time_evolution(n, t0+(m-1)*τg, τg, nsteps, rho0, H, G)
        append!(t_collected, t)
        append!(rho_collected, rho)
        rho0 = last(rho)  
    end

    rho0 = reset_qubit(rho0, rho_reset)

    t, rho = cdd_time_evolution(n, t0+M*τg, τg, nsteps, rho0, H, G)

    append!(t_collected, t)
    append!(rho_collected, rho)

    return t_collected, rho_collected
end


"""
    periodic_control_unconditional(M_max, n, t0, τg, nsteps, rho0, rho_reset, H, G, fout=nothing)

Unconditional periodic control implementation that averages over all possible measurement outcomes.
Returns either:
- If fout=nothing: Array of M_max states
- If fout is provided: Array of M_max fout values
"""
function periodic_control_unconditional(
    M_max::Int64,
    n::Int64,
    t0::Float64,
    τg::Float64,
    nsteps::Int64, 
    rho0::Operator,
    rho_reset::Operator,
    H::TimeDependentSum,
    G::Union{SVector,Vector},
    fout::Union{Function,Nothing}=nothing
    )

    final_states = Vector{Operator}(undef, M_max+1)
    current_state = rho0
    
    for m in 0:M_max
        final_states[m+1] = current_state
        current_state = cdd_final_state(n, t0+m*τg, τg, nsteps, current_state, H, G)
        current_state = reset_qubit(current_state, rho_reset)
    end

    if fout === nothing
        return final_states
    else
        return [fout(state) for state in final_states]
    end
end


"""
    periodic_control_conditional(M_max, n, t0, τg, nsteps, rho0, rho_reset, H, G, measurement_basis, fout=nothing)

Conditional periodic control that performs measurements and state updates based on measurement outcomes.
"""
function periodic_control_conditional(
    M_max::Int64,
    n::Int64,
    t0::Float64,
    τg::Float64,
    nsteps::Int64,
    rho0::Operator,
    rho_reset::Operator,
    H::TimeDependentSum,
    G::Union{SVector,Vector},
    measurement_basis::Operator,
    fout::Union{Function,Nothing}=nothing
    )
    
    # Get system dimensions to construct projectors
    b_sys = basis(measurement_basis)
    Id_sys = one(b_sys)
    
    # If there's a bath, get its dimension and identity
    has_bath = size(rho0)[1] > size(measurement_basis)[1]
    if has_bath
        b_bath = tensor([SpinBasis(1//2) for _ in 1:Int(log2(size(rho0)[1])/2-1)]...)
        Id_bath = one(b_bath)
    end
    
    # Store states and outcomes for this realization
    states = Vector{Operator}(undef, M_max+1)
    outcomes = Vector{Int64}(undef, M_max)
    
    # Initial state
    current_state = rho0
    states[1] = current_state
    
    for m in 1:M_max
        # Evolve state through CDD sequence
        current_state = cdd_final_state(n, t0+(m-1)*τg, τg, nsteps, current_state, H, G)
        
        # Calculate expectation value in measurement basis
        if has_bath
            rho_sys = ptrace(current_state, 2:Int(log2(size(current_state)[1])/2))
        else
            rho_sys = current_state
        end
        
        # Generate measurement outcome (binary: 0 or 1)
        expect_val = real(tr(measurement_basis * rho_sys))
        prob_zero = 0.5 * (1 + expect_val)
        outcome = rand() < prob_zero ? 0 : 1
        outcomes[m] = outcome
        
        # Update state based on measurement
        if has_bath
            proj_sys = 0.5 * (Id_sys + (outcome == 0 ? 1 : -1) * measurement_basis)
            projector = proj_sys ⊗ Id_bath
            current_state = projector * current_state * dagger(projector)
            normalize!(current_state)
            rho_bath = ptrace(current_state, 1)
            current_state = rho_reset ⊗ rho_bath
        else
            current_state = rho_reset
        end
        states[m+1] = current_state
    end
    
    if fout === nothing
        return states, outcomes
    else
        return [fout(state) for state in states], outcomes
    end
end


"""
    periodic_control(M_max, n, t0, τg, nsteps, rho0, rho_reset, H, G, measurement_scheme)

Main periodic control function that dispatches to either conditional or unconditional evolution
based on the measurement scheme configuration.

# Arguments
- `M_max::Int64`: The maximum number of control cycles
- `n::Int64`: The number of qubits in the system
- `t0::Float64`: The initial time
- `τg::Float64`: The duration of each control cycle
- `nsteps::Int64`: The number of time steps in each control cycle
- `rho0::Operator`: The initial density operator of the system (and bath if present)
- `rho_reset::Operator`: The reset state for the system part after each measurement
- `H::TimeDependentSum`: The time-dependent Hamiltonian of the system
- `G::Union{SVector,Vector}`: Fundamental control pulses
- `measurement_scheme::MeasurementScheme`: Configuration for measurements including:
    - measurement_basis: The operator defining the measurement basis (in system space)
    - conditional_evolution: Whether to perform conditional evolution
    - fout: Optional function to process output states

Bath handling:
- If a bath is present (rho0 dimension > rho_reset dimension):
    - Measurements are performed on the system part only
    - The system is reset to rho_reset while preserving bath state
    - Control pulses are automatically extended to the bath space

# Returns
- For conditional evolution (measurement_scheme.conditional_evolution = true):
    - Tuple of (states/fout_values, measurement_outcomes) where:
        - states is an array of M_max+1 states (if fout=nothing)
        - fout_values is an array of M_max+1 fout values (if fout is provided)
        - measurement_outcomes is an array of M_max binary values (0 or 1)
- For unconditional evolution:
    - Array of M_max+1 states (if fout=nothing)
    - Array of M_max+1 fout values (if fout is provided)
"""
function periodic_control(
    M_max::Int64,
    n::Int64,
    t0::Float64,
    τg::Float64,
    nsteps::Int64, 
    rho0::Operator,
    rho_reset::Operator,
    H::TimeDependentSum,
    G::Union{SVector,Vector},
    measurement_scheme::MeasurementScheme
    )
    
    if measurement_scheme.conditional_evolution
        if measurement_scheme.measurement_basis === nothing
            error("Measurement basis needed for conditional evolution")
        end
        return periodic_control_conditional(
            M_max, n, t0, τg, nsteps, rho0, rho_reset, H, G, 
            measurement_scheme.measurement_basis, measurement_scheme.fout
        )
    else
        return periodic_control_unconditional(M_max, n, t0, τg, nsteps, rho0, rho_reset, H, G, measurement_scheme.fout)
    end
end

end # module ControlledEvolution