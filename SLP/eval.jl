#= 
    @author: Thomas Chen
    @advisors: Yap Chee, Kai Hormann, Bingwei Zhang
    @date: 2025-06-18
    @description: Evaluates the given range of the slp under given assignment.
        The main function of this file is 
        evaluate_slp_range(slp, range_name, assignments)
            evaluates the given range of the slp under given assignment.

=# 

include("SLP.jl")


#=
    compute_operation(op::Symbol, left_val::Union{Float64, myInterval}, right_val::Union{Float64, myInterval})::Union{Float64, myInterval}
    
    Performs arithmetic operations on two operands based on the given operator symbol.
    
    # Arguments
    - `op::Symbol`: The operation symbol (:+, :-, :*, :^)
    - `left_val::Union{Float64, myInterval}`: Left operand value (can be Float64, myInterval, etc.)
    - `right_val::Union{Float64, myInterval}`: Right operand value (can be Float64, myInterval, etc.)
    
    # Returns
    - Union{Float64, myInterval}: The result of the arithmetic operation
    
    # Notes
    - For power operations (:^), the exponent must be a constant integer
    - If right_val is a myInterval for power operations, it must have equal lower and upper bounds(Meaning it represents an exact number)
    - Throws an error for unsupported operations or non-integer exponents in power operations
=#
function compute_operation(op::Symbol, left_val::Union{Float64, myInterval}, right_val::Union{Float64, myInterval})::Union{Float64, myInterval}
    if op == :+
        return left_val + right_val
    elseif op == :-
        return left_val - right_val
    elseif op == :*
        return left_val * right_val
    elseif op == :^
        if right_val isa myInterval
            if right_val.lower != right_val.upper
                error("Power operation requires constant exponent")
            end
            exponent = Int(right_val.lower)
        elseif right_val isa Float64 && isinteger(right_val)
            exponent = Int(right_val)
        else
            error("Power operation requires integer exponent")
        end
        left_val ^ exponent
    else
        error("Unsupported operation: $op")
    end
end


#=
    get_operand_value(slp::SLP, state::EvalState, idx::Union{Int, String})::Union{Float64, myInterval}
    
    Retrieves the value of an operand from an SLP using a separate evaluation workspace.
    
    # Arguments
    - `slp::SLP`: The structural SLP
    - `state::EvalState`: Per-evaluation runtime state
    - `idx::Union{Int, String}`: The index of the operand to retrieve
    
    # Returns
    - Union{Float64, myInterval}: The value associated with the given index or constant
    
    # Notes
    - Constants and implicit multiplication start with @ symbol
    - Negative indices refer to variables in slp.vars / state.var_values
    - Positive indices refer to intermediate results in slp.codelist / state.inst_values
    - Instruction values are cached in `state.inst_values` for the current evaluation only
    - Throws an error if the operand index is not found in the SLP
=#
function get_operand_value(slp::SLP, state::EvalState, idx::Union{Int, String})::Union{Float64, myInterval}
    if idx isa String
        if startswith(idx, "@")
            if contains(idx, "*") 
                parts = split(idx[2:end], "*")
                if length(parts) == 2
                    coeff = parse(Float64, parts[1])
                    var_idx = parse(Int, parts[2])
                    var_value = get_operand_value(slp, state, var_idx)
                    return coeff * var_value
                end
            else
                return parse(Float64, idx[2:end])
            end
        end
    elseif idx isa Int && idx < 0
        var_pos = -idx
        if var_pos <= length(state.var_values)
            return state.var_values[var_pos]
        end
    else
        if idx <= length(slp.codelist)
            cached_value = state.inst_values[idx]
            if cached_value !== nothing
                return cached_value
            else
                _, op, left_idx, right_idx = slp.codelist[idx]
                left_val = get_operand_value(slp, state, left_idx)
                right_val = get_operand_value(slp, state, right_idx)
                new_val = compute_operation(op, left_val, right_val)
                state.inst_values[idx] = new_val
                return new_val
            end
        end
    end
    error("Operand $idx not found in SLP")
end


#=
    load_assignments!(slp::SLP, state::EvalState, assignments::Dict{Symbol, T}) where T <: Union{Float64, myInterval}

    Loads variable assignments into the evaluation workspace and clears prior instruction values.
=#
function load_assignments!(slp::SLP, state::EvalState, assignments::Dict{Symbol, T}) where T <: Union{Float64, myInterval}
    if length(state.var_values) != length(slp.vars)
        error("EvalState variable buffer size does not match SLP variables")
    end
    if length(state.inst_values) != length(slp.codelist)
        error("EvalState instruction buffer size does not match SLP codelist")
    end

    for i in 1:length(slp.vars)
        _, variable = slp.vars[i]
        if haskey(assignments, variable)
            state.var_values[i] = assignments[variable]
        else
            error("Variable $variable not defined in assignmentsDict")
        end
    end

    reset!(state)
end


#=
    evaluate_slp_range(slp::SLP, state::EvalState, range_name::String, assignments::Dict{Symbol, T})::T where T <: Union{Float64, myInterval}
    
    Evaluates an SLP range (original function or derivative) with given variable values
    using a reusable evaluation workspace.
=#
function evaluate_slp_range(slp::SLP, state::EvalState, range_name::String, assignments::Dict{Symbol, T})::T where T <: Union{Float64, myInterval}
    if !haskey(slp.slp_ranges, range_name)
        error("Range $range_name does not exist in SLP")
    end

    range_start, range_end = slp.slp_ranges[range_name]

    if range_start > range_end 
        return T <: myInterval ? myInterval(0.0, 0.0) : 0.0
    end

    load_assignments!(slp, state, assignments)

    for i in range_start:range_end
        if state.inst_values[i] === nothing
            _, op, left_idx, right_idx = slp.codelist[i]
            left_val = get_operand_value(slp, state, left_idx)
            right_val = get_operand_value(slp, state, right_idx)
            state.inst_values[i] = compute_operation(op, left_val, right_val)
        end
    end

    result = state.inst_values[range_end]
    if result === nothing
        error("Failed to evaluate range $range_name")
    end
    return result
end


#=
    evaluate_slp_range(slp::SLP, range_name::String, assignments::Dict{Symbol, T})::T where T <: Union{Float64, myInterval}

    Convenience wrapper that allocates a fresh EvalState for a single evaluation.
=#
function evaluate_slp_range(slp::SLP, range_name::String, assignments::Dict{Symbol, T})::T where T <: Union{Float64, myInterval}
    state = EvalState(slp)
    return evaluate_slp_range(slp, state, range_name, assignments)
end
