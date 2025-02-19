module OutputProcessing

using QuantumOptics

export process_fidelity, process_expectation, process_fidelity_with_plus

"""
Calculate fidelity between states and a target state, handling bath tracing if needed.
Can process either a single state or an array of states against a fixed target.
"""
function process_fidelity(ρsb::Union{Operator, Vector{<:Operator}}, ρ_target::Operator)
    if ρsb isa Operator
        return _process_single_fidelity(ρsb, ρ_target)
    else
        return [_process_single_fidelity(state, ρ_target) for state in ρsb]
    end
end

"""
Internal helper for calculating fidelity of a single state
"""
function _process_single_fidelity(ρsb::Operator, ρ_target::Operator)
    if size(ρsb) == size(ρ_target)
        return real(tr(ρsb * ρ_target))
    end
    Nb = Int(log2(size(ρsb)[1])-1)
    ρs = ptrace(ρsb, 2:(Nb+1))
    return real(tr(ρs * ρ_target))
end

"""
Calculate expectation value of an operator for given states
"""
function process_expectation(operator::Operator, states)
    return real.(expect.(Ref(operator), states))
end

"""
Calculate fidelity with respect to the |+⟩ state for a list of states
"""
function process_fidelity_with_plus(states)
    b = SpinBasis(1//2)
    plus_state = normalize!(spinup(b) + spindown(b))
    target = projector(plus_state)
    return process_fidelity(states, target)
end

end # module OutputProcessing