#= 
    @author: Thomas Chen
    @advisors: Yap Chee, Kai Hormann, Bingwei Zhang
    @date: 2025-06-18
    @description: Evaluates the given range of the slp under given assignment.
        The main function of this file is 
        evaluate_slp_range(slp, range_name, assignments)
            evaluates the given range of the slp under given assignment.

=# 


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
    get_operand_value(slp::SLP, idx::Union{Int, String})::Union{Float64, myInterval}
    
    Retrieves the value of an operand from an SLP (Straight Line Program) structure.
    
    # Arguments
    - `slp::SLP`: The SLP structure containing variables, constants, and code list
    - `idx::Union{Int, String}`: The index of the operand to retrieve
    
    # Returns
    - Union{Float64, myInterval}: The value associated with the given index, or constant
    
    # Notes
    - Constants and implicit multiplication starts with @ symbol
    - Negative indices refer to variables in slp.vars
    - Positive indices refer to intermediate results in slp.codelist
    - If a value is not precomputed (equals myInterval(0.0, 0.0)), it recursively
      computes the value using the operation and operand indices
    - Throws an error if the operand index is not found in the SLP
=#
function get_operand_value(slp::SLP, idx::Union{Int, String})::Union{Float64, myInterval}
    if idx isa String
        if startswith(idx, "@")
            if contains(idx, "*") 
                parts = split(idx[2:end], "*")
                if length(parts) == 2
                    coeff = parse(Float64, parts[1])
                    var_idx = parse(Int, parts[2])
                    var_value = get_operand_value(slp, var_idx)
                    return coeff * var_value
                end
            else
                return parse(Float64, idx[2:end])
            end
        end
    elseif idx isa Int && idx < 0 # Variable, get its value
        return slp.vars[-idx][3]
    else
        if idx <= length(slp.codelist)
            _, op, left_idx, right_idx, value = slp.codelist[idx]
            if value !== myInterval(0.0, 0.0)
                return value
            else
                left_val = get_operand_value(slp, left_idx)
                right_val = get_operand_value(slp, right_idx)
                return compute_operation(op, left_val, right_val)
            end
        end
    end
    error("Operand $idx not found in SLP")
end


#=
    evaluate_slp_range(slp::SLP, range_name::String, assignments::Dict{Symbol, T})::T where T <: Union{Float64, myInterval}
    
    Evaluates an SLP range (original function or derivative) with given variable values.
    
    Arguments:
    - slp: The SLP structure containing the computational graph
    - range_name: Name of the range to evaluate (e.g., "", "x", "x^2", "xy", etc.)
    - assignments: Dictionary mapping variable symbols to their values (Float64 or myInterval)
    
    Returns:
    - The evaluated result as the same type as the input values
=#
function evaluate_slp_range(slp::SLP, range_name::String, assignments::Dict{Symbol, T})::T where T <: Union{Float64, myInterval}
    if !haskey(slp.slp_ranges, range_name)
        error("Range $range_name does not exist in SLP")
    end

    range_start, range_end = slp.slp_ranges[range_name]

    # range_end = min(range_end, length(slp.codelist))
    # range_start = max(range_start, 1)

    if range_start > range_end 
        return T <: myInterval ? myInterval(0.0, 0.0) : 0.0
    end
    
    for i in 1:length(slp.vars)
        index, variable, _ = slp.vars[i]
        if variable isa Symbol && haskey(assignments, variable)
            computed_value = assignments[variable]
            slp.vars[i] = (index, variable, computed_value)
        else
            error("Variable $variable not defined in assignmentsDict")
        end
    end

    for i in range_start:range_end
        inst_idx, op, left_idx, right_idx, _ = slp.codelist[i]
        left_val = get_operand_value(slp, left_idx)
        right_val = get_operand_value(slp, right_idx)
        new_val = compute_operation(op, left_val, right_val)
        slp.codelist[i] = (inst_idx, op, left_idx, right_idx, new_val)
    end

    return slp.codelist[range_end][end]
end
