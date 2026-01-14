#= 
    @author: Thomas Chen
    @advisors: Chee Yap, Kai Hormann, Bingwei Zhang
    @date: 2025-06-17
    @description: SLP struct
		The SLP (Straight Line Program) struct represents a computational
		graph (see details below).
            The useful functions in this file are:
			-- get_constant_value(idx) - Returns the constant value as
				represented by a string.
                - E.g. "@2.5" -> 2.5 
			-- is_constant_value(idx) - Returns true if the operand is a
				constant value.
			-- get_or_add_variable(slp, var) - Returns the index of a
				variable in SLP. If the variable does not exist, it is
				added.
			-- create_variable_instruction!(slp, var) - Creates a
				variable instruction in SLP.
		MISSING:
			What are the primitive operations in the codelists??
=# 

# Type alias for numerical values throughout the SLP structure
# - Float64: standard floating-point numbers for exact computations
# - myInterval: custom interval type for interval arithmetic computations
const valueType = Union{Float64, myInterval}

#= 
The SLP (Straight Line Program) struct represents a computational graph
	that stores mathematical operations and their dependencies in a
	linear format.  It supports symbolic differentiation by maintaining
	shared storage for original functions and their derivatives.
    
    Fields:
    - vars: Vector storing variables with negative indices (-1 to -∞)
                 Each tuple contains (index, variable, value)
				 					?? variable name like x, y?? 
	- codelist: Vector storing computational operations with positive
		indices (1 to ∞) Each tuple contains (index, operation_symbol,
		operand1_idx, operand2_idx, result_value)
	- slp_ranges: Dictionary mapping SLP names to their index ranges in
								?? SLP name ??
		codelist Enables tracking of original SLP and derivative SLPs
		within shared storage
			Convention: "" represents the original function f,
				 "x^ky^l" represents the derivative of f w.r.t to x k
				 times and y l times. 
	- global_deriv_map: Dictionary mapping (instruction_index, variable)
	  	pairs to their derivative instruction index
			  Supports efficient retrieval of derivative instructions for
			  given variables
=#
struct SLP 
    vars::Vector{Tuple{Int, Symbol, valueType}}
    codelist::Vector{Tuple{Int, Symbol, Union{Int, String},
						Union{Int, String}, valueType}}
    slp_ranges::Dict{String, Tuple{Int, Int}}
    global_deriv_map::Dict{Tuple{Int, Symbol}, Union{Int, String}}  
end


#=
    get_constant_value(idx::Union{Int, String})::Union{Float64, Nothing}
    
	Extracts constant value from SLP operand. Returns nothing if operand
		is not a constant.
    
    # Arguments
    - `slp::SLP`: SLP struct instance
	- `idx::Union{Int, String}`: Index or string representation of the
	  	operand
    
    # Returns
	- `Union{Float64, Nothing}`: Constant value if operand is a constant,
	  nothing otherwise
=#
function get_constant_value(idx::Union{Int, String})::Union{Float64, Nothing}
    if idx isa String && startswith(idx, "@")
        if contains(idx, "*")
            return nothing
        else
            return parse(Float64, idx[2:end])
        end
    end
    return nothing
end

#=
    is_constant_value(idx::Union{Int, String}, value::Float64)::Bool
    
	Checks if operand is a specific constant value. Returns false if
			operand is not a constant.
    
    # Arguments
    - `slp::SLP`: SLP struct instance
	- `idx::Union{Int, String}`: Index or string representation of the
	  		operand
    - `value::Float64`: Constant value to compare
    
    # Returns
	- `Bool`: True if operand is a constant and matches the value, false
	  		otherwise
=#
function is_constant_value(idx::Union{Int, String}, value::Float64)::Bool
    if idx isa String && startswith(idx, "@")
        if contains(idx, "*")
            return false
        else
            const_val = parse(Float64, idx[2:end])
            return const_val == value
        end
    end
    return false
end

#=
    get_or_add_variable(slp::SLP, var::Symbol)::Int
    
    Helper function to get existing variable or add new one to SLP
    
    # Arguments
    - `slp::SLP`: SLP struct instance
    - `var::Symbol`: Variable name to be added or retrieved
    
    # Returns
    - `Int`: Index of the variable in SLP
=#
function get_or_add_variable(slp::SLP, var::Symbol)::Int
    for (idx, val, _) in slp.vars 
        if val isa Symbol && val == var
            return idx
        end
    end
    
    new_idx = isempty(slp.vars) ? -1 : minimum([idx for (idx, _, _) in slp.vars]) - 1
    push!(slp.vars, (new_idx, var, myInterval(0.0, 0.0)))
    return new_idx
end

#= 
    create_variable_instruction!(slp::SLP, var::Symbol)::Int
    Helper function to create a variable instruction in the SLP
    
    # Arguments
    - `slp::SLP`: SLP struct instance
    - `var::Symbol`: Variable name to be added or retrieved
    
    # Returns
    - `Int`: Index of the variable in SLP
=#
function create_variable_instruction!(slp::SLP, var::Symbol)::Int
    var_idx = get_or_add_variable(slp, var)
    zero_literal = "@0.0" 
    
    for (i, instruction) in enumerate(slp.codelist)
        out_idx, op, left_idx, right_idx, _ = instruction
        if op == :+ && ((left_idx == var_idx
					&& right_idx == zero_literal) || 
                       (left_idx == zero_literal && right_idx == var_idx))
            return out_idx
        end
    end
    
    new_out_idx = length(slp.codelist) + 1
    push!(slp.codelist,
		(new_out_idx, :+, var_idx, zero_literal, myInterval(0.0, 0.0)))
    return new_out_idx
end
