module ControlledEvolution

using QuantumOpticsBase, QuantumOptics, StaticArrays, LinearAlgebra, Distributions, Random
using ..CoreTypes

export cdd, cdd_filter_function, reset_qubit, cdd_final_state, cdd_time_evolution, periodic_control_final_state, 
periodic_control_time_evolution, gate_fidelity, periodic_control, periodic_control_conditional


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


"""
periodic_control(M_max, n, t0, τg, nsteps, rho0, rho_reset, H, G)

This function calculates the final states of a quantum system under periodic control.

# Arguments
- `M_max::Int64`: The maximum number of control cycles.
- `n::Int64`: The number of qubits in the system.
- `t0::Float64`: The initial time.
- `τg::Float64`: The duration of each control cycle.
- `nsteps::Int64`: The number of time steps in each control cycle.
- `rho0::Operator`: The initial density operator of the system.
- `rho_reset::Operator`: The reset density operator after each control cycle.
- `H::TimeDependentSum`: The time-dependent Hamiltonian of the system.
- `G::Union{SVector,Vector}`: Fundamental control pulses.

# Returns
- `final_states`: An array of final density operators after each control cycle
    and idling gate implementation.

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
    G::Union{SVector,Vector}
    )

    # final states after m repetitions of the control sequence
    ctrl_states = Vector{Operator}(undef, M_max+1)
    ctrl_states[1] = rho0
    for m in 1:M_max
        rho0 = cdd_final_state(n, t0+(m-1)*τg, τg, nsteps, rho0, H, G)
        ctrl_states[m+1] = rho0
    end

    final_states = Vector{Operator}(undef, M_max+1)
    for m in 0:M_max
        reset_state = reset_qubit(ctrl_states[m+1], rho_reset)
        final_state = cdd_final_state(n, t0+m*τg, τg, nsteps, reset_state, H, G)
        final_states[m+1] =  final_state
    end

    return final_states

end


# returns the final state after the whole protocol
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


function gate_fidelity(ρsb, ρ_target)
    # If dimensions match, there's no bath to trace out
    if size(ρsb) == size(ρ_target)
        return real(tr(ρsb * ρ_target))
    end
    Nb = Int(log2(size(ρsb)[1])-1)
    ρs = ptrace(ρsb, 2:(Nb+1))
    return real(tr(ρs * ρ_target))
end


"""
    periodic_control_conditional(M_max, n, t0, τg, nsteps, rho0, rho_reset, H, G, measurement_basis)

This function implements periodic control with conditional evolution based on measurement outcomes.

# Arguments
- `M_max::Int64`: The maximum number of control cycles
- `n::Int64`: The number of qubits in the system
- `t0::Float64`: The initial time
- `τg::Float64`: The duration of each control cycle
- `nsteps::Int64`: The number of time steps in each control cycle
- `rho0::Operator`: The initial density operator of the system
- `rho_reset::Operator`: The reset state for the system part after each measurement
- `H::TimeDependentSum`: The time-dependent Hamiltonian of the system
- `G::Union{SVector,Vector}`: Fundamental control pulses
- `measurement_basis::Operator`: The operator defining the measurement basis (in system space)

# Returns
- `expectation_values`: Array of expectation values ⟨σⱼ⟩ after each cycle
- `measurement_outcomes`: Array of measurement outcomes (0 or 1) for each cycle
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
    measurement_basis::Operator
    )
    
    # Initialize arrays to store results
    expectation_values = Float64[]
    measurement_outcomes = Int64[]
    
    # Current state
    current_state = rho0
    
    # Get system dimensions to construct projectors
    b_sys = basis(measurement_basis)
    Id_sys = one(b_sys)
    
    # If there's a bath, get its dimension and identity
    has_bath = size(current_state)[1] > size(measurement_basis)[1]
    if has_bath
        b_bath = tensor([SpinBasis(1//2) for _ in 1:Int(log2(size(current_state)[1])/2-1)]...)
        Id_bath = one(b_bath)
    end
    
    for m in 1:M_max
        # Evolve state through CDD sequence
        current_state = cdd_final_state(n, t0+(m-1)*τg, τg, nsteps, current_state, H, G)
        
        # Calculate expectation value in measurement basis
        if has_bath
            # Partial trace over bath to get system state for expectation value
            rho_sys = ptrace(current_state, 2:Int(log2(size(current_state)[1])/2))
        else
            rho_sys = current_state
        end
        
        expect_val = real(tr(measurement_basis * rho_sys))
        push!(expectation_values, expect_val)
        
        # Calculate measurement probability (P(0) = 1/2(1 + ⟨σⱼ⟩))
        prob_zero = 0.5 * (1 + expect_val)
        
        # Generate random measurement outcome
        outcome = rand() < prob_zero ? 0 : 1
        push!(measurement_outcomes, outcome)
        
        if has_bath
            # Construct appropriate projector based on outcome
            proj_sys = 0.5 * (Id_sys + (outcome == 0 ? 1 : -1) * measurement_basis)
            projector = proj_sys ⊗ Id_bath
            
            # Apply projection and normalize
            current_state = projector * current_state * dagger(projector)
            normalize!(current_state)
            
            # Get bath state after measurement
            rho_bath = ptrace(current_state, 1)
            
            # Prepare next state as tensor product of reset system state and current bath state
            current_state = rho_reset ⊗ rho_bath
        else
            # For no bath case, just reset the state
            current_state = rho_reset
        end
    end
    
    return expectation_values, measurement_outcomes
end


end # module ControlledEvolution