#= 
    @author: Thomas Chen
    @advisors: Yap Chee, Kai Hormann, Bingwei Zhang
    @date: 2025-06-18
    @description: Derivative Computation for SLP polynomials
        Two main functions: 
            - compute_all_derivatives!(poly::Polynomial)
			- ???
		What are the primitive operations in the codelists??
=# 

include("polynomial.jl")
include("SLP.jl")
include("myInterval.jl")
include("utils.jl")

COEFF_PARSE_CACHE = Dict{String, Tuple{Float64, Int}}()


#=
	get_or_compute_derivative(poly::Polynomial, idx::Union{Int, String},
			var::Symbol)::Union{Int, String}
    
    Gets the derivative of instruction at idx with respect to var.
    If not already computed, recursively computes it.
    
    # Arguments
    - `poly::Polynomial`: The polynomial whose derivative is being computed
	- `idx::Union{Int, String}`: The index of the instruction in the
	  polynomial's SLP, or a constant value
	- `var::Symbol`: The variable with respect to which the derivative is
	  being computed
    
    # Returns
	- `::Union{Int, String}`: The index of the derivative instruction in
	  the polynomial's SLP or a string representing the the constant
	  derivative
    
    # Notes
	- If the derivative has already been computed, it is returned from
	  the global derivative map
    - If the instruction is a constant, the derivative is 0
	- If the instruction is a variable, the derivative is 1 if the
	  variable matches, and 0 otherwise
	- If the instruction is an operation, the derivative is computed
	  based on the operation and its operands
=#
function get_or_compute_derivative(poly::Polynomial,
		idx::Union{Int, String}, var::Symbol)::Union{Int, String}
    global COEFF_PARSE_CACHE
    if idx isa String
        if startswith(idx, "@")
            if contains(idx, "*")
                # Use cached parsing
                coeff, var_idx = if haskey(COEFF_PARSE_CACHE, idx)
                    COEFF_PARSE_CACHE[idx]
            else
                    content = idx[2:end] 
                    parts = split(content, "*")
                    if length(parts) == 2
                        coeff_val = parse(Float64, parts[1])
                        var_idx_val = parse(Int, parts[2])
                        COEFF_PARSE_CACHE[idx] = (coeff_val, var_idx_val)
                        (coeff_val, var_idx_val)
                    else
                        error("Invalid coefficient string format: $idx")
                    end
            end
                
            var_deriv = get_or_compute_derivative(poly, var_idx, var)
                
                if var_deriv == "@0.0"
                    return "@0.0"
                end
                
                if var_deriv == "@1.0"
                    return "@$(coeff)"
                end
                
                return add_optimized_instruction!(poly.slp, :*,
					"@$(coeff)", var_deriv)
            else
                return "@0.0"
            end
        end
    end
    
    if haskey(poly.slp.global_deriv_map, (idx, var))
        return poly.slp.global_deriv_map[(idx, var)]
    end
    
    if idx < 0
        try
            _, var_symbol, _ = poly.slp.vars[abs(idx)]
            if var_symbol == var
                poly.slp.global_deriv_map[(idx, var)] = "@1.0"
                return "@1.0"
            else
                poly.slp.global_deriv_map[(idx, var)] = "@0.0"
                return "@0.0"
            end
        catch BoundsError
            error("Variable with index $idx not found in SLP")
        end
    end

    if idx > 0 && idx <= length(poly.slp.codelist)
        instruction = poly.slp.codelist[idx]
        out_idx, op, left_idx, right_idx = instruction
        
        if op == :+
            left_deriv_idx = get_or_compute_derivative(poly, left_idx, var)
            right_deriv_idx = get_or_compute_derivative(poly, right_idx, var)
            result_idx = add_optimized_instruction!(poly.slp,
				:+, left_deriv_idx, right_deriv_idx)
            poly.slp.global_deriv_map[(idx, var)] = result_idx
            return result_idx
            
        elseif op == :-
            left_deriv_idx = get_or_compute_derivative(poly, left_idx, var)
            right_deriv_idx = get_or_compute_derivative(poly, right_idx, var)
            result_idx = add_optimized_instruction!(poly.slp,
				:-, left_deriv_idx, right_deriv_idx)
            poly.slp.global_deriv_map[(idx, var)] = result_idx
            return result_idx
            
        elseif op == :*
            left_deriv_idx = get_or_compute_derivative(poly, left_idx, var)
            right_deriv_idx = get_or_compute_derivative(poly, right_idx, var)
            
            term1_idx = add_optimized_instruction!(poly.slp, :*, left_deriv_idx, right_idx)
            term2_idx = add_optimized_instruction!(poly.slp, :*, left_idx, right_deriv_idx)
            result_idx = add_optimized_instruction!(poly.slp, :+, term1_idx, term2_idx)
            poly.slp.global_deriv_map[(idx, var)] = result_idx
            return result_idx
            
        elseif op == :^
            exp_val = get_constant_value(right_idx)
            if exp_val === nothing
                error("Cannot differentiate non-constant exponent")
            end
            
            if exp_val == 0.0
                poly.slp.global_deriv_map[(idx, var)] = "@0.0"
                return "@0.0"
            elseif exp_val == 1.0
                base_deriv_idx = get_or_compute_derivative(poly, left_idx, var)
                poly.slp.global_deriv_map[(idx, var)] = base_deriv_idx
                return base_deriv_idx
            end

            base_deriv_idx = get_or_compute_derivative(poly, left_idx, var)
            if is_constant_value(base_deriv_idx, 0.0)
                poly.slp.global_deriv_map[(idx, var)] = "@0.0"
                return "@0.0"  
            end

            n_idx = "@$(exp_val)"
            n_minus_1_idx = "@$(exp_val - 1.0)"
            
            power_idx = add_optimized_instruction!(poly.slp, :^, left_idx, n_minus_1_idx)
            coeff_idx = add_optimized_instruction!(poly.slp, :*, n_idx, power_idx)
            result_idx = add_optimized_instruction!(poly.slp, :*, coeff_idx, base_deriv_idx)
            poly.slp.global_deriv_map[(idx, var)] = result_idx
            return result_idx
        else
            error("Unsupported operation: $op")
        end
    end
    
    poly.slp.global_deriv_map[(idx, var)] = "@0.0"
    return "@0.0"
end


#=
    create_derivative_mapping(slp::SLP, diff_var::Symbol)
    
    Creates a mapping from original instruction indices to their derivative indices for a given variable.
    
    # Arguments
    - `slp::SLP`: The SLP structure containing variables, constants, and code list
    - `diff_var::Symbol`: The variable with respect to which the derivative is being computed
    
    # Notes
    - If the derivative has already been computed, it is returned from the global derivative map
    - If the instruction is a constant, the derivative is 0
    - If the instruction is a variable, the derivative is 1 if the variable matches, and 0 otherwise
    - If the instruction is an operation, the derivative is computed based on the operation and its operands
=#
function create_derivative_mapping(slp::SLP, diff_var::Symbol)
    if slp.global_deriv_map === nothing
        slp.global_deriv_map = Dict{Tuple{Int, Symbol}, Union{Int, String}}()
    end
    
    for (idx, val) in slp.vars
        key = (idx, diff_var)
        if !haskey(slp.global_deriv_map, key)
            if val isa Symbol && val == diff_var
                slp.global_deriv_map[key] = "@1.0"
            else
                slp.global_deriv_map[key] = "@0.0"
            end
        end
    end
end

#=
    add_optimized_instruction!(slp::SLP, op::Symbol, left_idx::Union{Int, String}, right_idx::Union{Int, String})::Union{Int, String}

    Helper function to optimize operations by skipping unnecessary instructions
        
    # Arguments
    - `slp::SLP`: The SLP structure containing variables, constants, and code list
    - `op::Symbol`: The operation symbol (+, -, *, ^)
    - `left_idx::Union{Int, String}`: The index of the left operand in the code list or a constant value
    - `right_idx::Union{Int, String}`: The index of the right operand in the code list or a constant value
    
    # Returns
    - `::Union{Int, String}`: The index of the optimized instruction in the code list or a constant value
    
    # Notes
    - If the operation is addition and the left operand is a constant 0, the right operand is returned
    - If the operation is subtraction and the right operand is a constant 0, the left operand is returned
    - If the operation is multiplication and either operand is a constant 0, the constant 0 is returned
    - If the operation is multiplication and either operand is a constant 1, the other operand is returned
    - If the operation is exponentiation and the right operand is a constant 1, the left operand is returned
    - If the operation is exponentiation and the right operand is a constant 0, the constant 1 is returned
=#
function add_optimized_instruction!(slp::SLP, op::Symbol, left_idx::Union{Int, String}, right_idx::Union{Int, String})::Union{Int, String}
    global INSTRUCTION_HASH
    
    if op == :+
        if is_constant_value(left_idx, 0.0)
            return right_idx
        elseif is_constant_value(right_idx, 0.0)
            return left_idx
        end
    
    elseif op == :-
        if is_constant_value(right_idx, 0.0)
            return left_idx
        end
    
    elseif op == :*
        if is_constant_value(left_idx, 0.0) || is_constant_value(right_idx, 0.0)
            return "@0.0"
        end
        if is_constant_value(left_idx, 1.0)
            return right_idx
        elseif is_constant_value(right_idx, 1.0)
            return left_idx
        end
        
        # Handle implicit multiplication generation for derivatives
        # Case 1: constant * variable/instruction
        if left_idx isa String && startswith(left_idx, "@") && !contains(left_idx, "*") && right_idx isa Int
            coeff = parse(Float64, left_idx[2:end])
            return "@$(coeff)*$(right_idx)"
        end
        
        # Case 2: variable/instruction * constant
        if right_idx isa String && startswith(right_idx, "@") && !contains(right_idx, "*") && left_idx isa Int
            coeff = parse(Float64, right_idx[2:end])
            return "@$(coeff)*$(left_idx)"
        end
        
        # Case 3: implicit multiplication * constant (e.g., @2.0*1 * @3.0 = @6.0*1)
        if left_idx isa String && startswith(left_idx, "@") && contains(left_idx, "*") && 
           right_idx isa String && startswith(right_idx, "@") && !contains(right_idx, "*")
            # Parse left implicit multiplication
            content = left_idx[2:end]
            parts = split(content, "*")
            if length(parts) == 2
                left_coeff = parse(Float64, parts[1])
                var_idx = parse(Int, parts[2])
                right_coeff = parse(Float64, right_idx[2:end])
                combined_coeff = left_coeff * right_coeff
                return "@$(combined_coeff)*$(var_idx)"
            end
        end
        
        # Case 4: constant * implicit multiplication (e.g., @3.0 * @2.0*1 = @6.0*1)
        if left_idx isa String && startswith(left_idx, "@") && !contains(left_idx, "*") &&
           right_idx isa String && startswith(right_idx, "@") && contains(right_idx, "*")
            # Parse right implicit multiplication
            content = right_idx[2:end]
            parts = split(content, "*")
            if length(parts) == 2
                right_coeff = parse(Float64, parts[1])
                var_idx = parse(Int, parts[2])
                left_coeff = parse(Float64, left_idx[2:end])
                combined_coeff = left_coeff * right_coeff
                return "@$(combined_coeff)*$(var_idx)"
            end
        end
    elseif op == :^
        if is_constant_value(right_idx, 1.0)
            return left_idx
        end
        if is_constant_value(right_idx, 0.0)
            return "@1.0"
        end
    end
    
    instruction_key = (op, left_idx, right_idx)
    if haskey(INSTRUCTION_HASH, instruction_key)
        return INSTRUCTION_HASH[instruction_key]
    end
    
    # Check commutative operations (+ and *)
    if op == :+ || op == :*
        commutative_key = (op, right_idx, left_idx)
        if haskey(INSTRUCTION_HASH, commutative_key)
            return INSTRUCTION_HASH[commutative_key]
        end
    end

    new_out_idx = length(slp.codelist) + 1
    push!(slp.codelist, (new_out_idx, op, left_idx, right_idx))

    INSTRUCTION_HASH[instruction_key] = new_out_idx
    if op == :+ || op == :*
        INSTRUCTION_HASH[(op, right_idx, left_idx)] = new_out_idx
    end

    return new_out_idx
end

#=
    compute_derivative_slp!(poly::Polynomial, var::Symbol, existing_deriv::Dict{Symbol, Int} = Dict{Symbol, Int}())::String

    Computes the derivative of the polynomial with respect to var and adds it to the SLP
    Assumes the base SLP (original function or previous derivative) is already computed

    # Arguments
    - `poly::Polynomial`: The polynomial structure containing the SLP and derivative information
    - `var::Symbol`: The variable with respect to which the derivative is being computed
	- `existing_deriv::Dict{Symbol, Int}`: A mapping of variables to
	  		their derivative orders
    
    # Returns
    - `::String`: The key for the new derivative in the SLP ranges dictionary
    
    # Notes
    - If the derivative has already been computed, the existing key is returned
    - If the base SLP is missing, an error is thrown
    - The derivative is computed by incrementing the derivative order for the given variable
=#
function compute_derivative_slp!(poly::Polynomial, var::Symbol,
		existing_deriv::Dict{Symbol, Int} = Dict{Symbol, Int}())::String
    if poly.slp === nothing
        error("Polynomial must have an SLP to compute derivatives")
    end

    new_deriv = copy(existing_deriv)
    new_deriv[var] = get(new_deriv, var, 0) + 1
    new_key = combine_derivative_strings(new_deriv)

    if haskey(poly.slp.slp_ranges, new_key)
        return new_key
    end

    base_key = combine_derivative_strings(existing_deriv)
    if !haskey(poly.slp.slp_ranges, base_key)
        error("Base derivative $base_key not found")
    end

    base_range = poly.slp.slp_ranges[base_key]
    start_idx = length(poly.slp.codelist) + 1
    last_instruction_idx = start_idx - 1  

    #If derivative of base is already 0, then current derivative also 0
    if base_range[2] < base_range[1]
        poly.slp.slp_ranges[new_key] = (start_idx, start_idx - 1)
        return new_key
    end

    if poly.slp.global_deriv_map === nothing
        poly.slp.global_deriv_map = Dict{Tuple{Int, Symbol},
				Union{Int, String}}()
    end

    create_derivative_mapping(poly.slp, var)

    for i in base_range[1]:base_range[2]
        if i > length(poly.slp.codelist)
            break
        end

        instruction = poly.slp.codelist[i]
        out_idx, op, left_idx, right_idx = instruction
        
        if op == :+ 
            left_deriv_idx = get_or_compute_derivative(poly, left_idx, var)
            right_deriv_idx = get_or_compute_derivative(poly, right_idx, var)
            
            result_idx = add_optimized_instruction!(poly.slp,
					:+, left_deriv_idx, right_deriv_idx)
            poly.slp.global_deriv_map[(out_idx, var)] = result_idx

            if result_idx isa Int && result_idx > last_instruction_idx
                last_instruction_idx = result_idx
            end
            
        elseif op == :- 
            left_deriv_idx = get_or_compute_derivative(poly, left_idx, var)
            right_deriv_idx = get_or_compute_derivative(poly, right_idx, var)
            
            result_idx = add_optimized_instruction!(poly.slp,
					:-, left_deriv_idx, right_deriv_idx)
            poly.slp.global_deriv_map[(out_idx, var)] = result_idx
            
            if result_idx isa Int && result_idx > last_instruction_idx
                last_instruction_idx = result_idx
            end
            
        elseif op == :*
            left_deriv_idx = get_or_compute_derivative(poly, left_idx, var)
            right_deriv_idx = get_or_compute_derivative(poly, right_idx, var)
            
            term1_idx = add_optimized_instruction!(poly.slp, :*, left_deriv_idx, right_idx)
            term2_idx = add_optimized_instruction!(poly.slp, :*, left_idx, right_deriv_idx)
            result_idx = add_optimized_instruction!(poly.slp, :+, term1_idx, term2_idx)
            poly.slp.global_deriv_map[(out_idx, var)] = result_idx
            
            if result_idx isa Int && result_idx > last_instruction_idx
                last_instruction_idx = result_idx
            end
            
        elseif op == :^
            base_deriv_idx = get_or_compute_derivative(poly, left_idx, var)
            
            exp_val = get_constant_value(right_idx)
            if exp_val === nothing
                error("Cannot differentiate non-constant exponent")
            end
            
            if exp_val == 0.0
                poly.slp.global_deriv_map[(out_idx, var)] = "@0.0"
            elseif exp_val == 1.0
                poly.slp.global_deriv_map[(out_idx, var)] = base_deriv_idx
            else
                if is_constant_value(base_deriv_idx, 0.0)
                    poly.slp.global_deriv_map[(out_idx, var)] = "@0.0"
                else
                    n_idx = "@$exp_val"
                    n_minus_1_idx = "@$(exp_val - 1.0)"
                    
                    power_idx = add_optimized_instruction!(poly.slp, :^, left_idx, n_minus_1_idx)
                    coeff_idx = add_optimized_instruction!(poly.slp, :*, n_idx, power_idx)  
                    final_idx = add_optimized_instruction!(poly.slp, :*, coeff_idx, base_deriv_idx)
                    poly.slp.global_deriv_map[(out_idx, var)] = final_idx
                    if final_idx isa Int
                        if final_idx > last_instruction_idx
                            last_instruction_idx = final_idx
                        end
                    end
                end
            end
        end 
    end


    if !isempty(poly.slp.global_deriv_map)
        final_out_idx = base_range[2]
        if haskey(poly.slp.global_deriv_map, (final_out_idx, var))
            final_result_idx = poly.slp.global_deriv_map[(final_out_idx, var)]

            const_val = get_constant_value(final_result_idx)
            if const_val !== nothing && const_val == 0.0
                poly.slp.slp_ranges[new_key] = (start_idx, start_idx - 1)
            elseif final_result_idx isa Int && final_result_idx > 0 && final_result_idx < start_idx
                poly.slp.slp_ranges[new_key] = (final_result_idx, final_result_idx)
            else
                if last_instruction_idx < start_idx
                    if const_val !== nothing
                        dummy_idx = length(poly.slp.codelist) + 1
                        push!(poly.slp.codelist, (dummy_idx, :+, "@$(const_val)", "@0.0"))
                        poly.slp.slp_ranges[new_key] = (dummy_idx, dummy_idx)
                    else
                        if final_result_idx isa String
                            dummy_idx = length(poly.slp.codelist) + 1
                            push!(poly.slp.codelist, (dummy_idx, :+, final_result_idx, "@0.0"))
                            poly.slp.slp_ranges[new_key] = (dummy_idx, dummy_idx)
                        else
                            for (idx, val) in poly.slp.vars
                                if idx == final_result_idx && val isa Symbol
                                    dummy_idx = create_variable_instruction!(poly.slp, val)
                                    poly.slp.slp_ranges[new_key] = (dummy_idx, dummy_idx)
                                    break
                                end
                            end
                        end
                    end
                else
                    poly.slp.slp_ranges[new_key] = (start_idx, last_instruction_idx)
                end
            end
        else
            poly.slp.slp_ranges[new_key] = (start_idx, start_idx - 1)
        end
    else
        poly.slp.slp_ranges[new_key] = (start_idx, start_idx - 1)
    end
    
    return new_key
end



#=  
    compute_all_derivatives!(poly::Polynomial)
    Function to compute all derivatives up to max_order for each variable
    
    # Arguments
	- `poly::Polynomial`: The polynomial structure containing the SLP and
	  			derivative information
    
    # Notes
    - If the polynomial does not have an SLP, an error is thrown
    - If the original function SLP range is missing, an error is thrown
    - The function computes derivatives for each variable up to max_order
=# 

function compute_all_derivatives!(poly::Polynomial)
    empty!(INSTRUCTION_HASH)
    if poly.slp === nothing
        error("Polynomial must have an SLP to compute derivatives")
    end

    if !haskey(poly.slp.slp_ranges, "")
        error("Original function SLP range not found. The base SLP must be computed first.")
    end

    variables = get_variables(poly)
    
    to_compute = Vector{Dict{Symbol, Int}}()
    computed = Set{String}()

    push!(computed, "")

    for var in variables
        if get_max_order(poly, var) >= 1
            first_deriv = Dict{Symbol, Int}(var => 1)
            push!(to_compute, first_deriv)
        end
    end

    while !isempty(to_compute) 		# main loop!!!
        current_deriv = popfirst!(to_compute)
        current_key = combine_derivative_strings(current_deriv)

        if current_key in computed
            continue
        end

        base_deriv = copy(current_deriv)
        diff_var = nothing

        for (var, order) in current_deriv
            if order > 0
                base_deriv[var] = order - 1
                if base_deriv[var] == 0
                    delete!(base_deriv, var)
                end
                diff_var = var
                break
            end
        end

        if diff_var === nothing
            continue
        end

        base_key = combine_derivative_strings(base_deriv)

        if base_key in computed
            try 
                compute_derivative_slp!(poly, diff_var, base_deriv)
                push!(computed, current_key)
                
                for var in variables 
                    new_deriv = copy(current_deriv)
                    current_order = get(new_deriv, var, 0)
                    max_order = get_max_order(poly, var)

                    if current_order < max_order
                        new_deriv[var] = current_order + 1
                        new_key = combine_derivative_strings(new_deriv)

                        if !(new_key in computed) && !(new_deriv in to_compute)
                            push!(to_compute, new_deriv)
                        end
                    end
                end

            catch e
                println("Error computing derivative for $current_key: $e")
            end
        else 
            push!(to_compute, current_deriv)
        end
    end
end

#= 
    find_valid_base_derivative(target_deriv::Dict{Symbol, Int}, computed::Set{String})::Tuple{Dict{Symbol, Int}, Symbol}
    Function to find the valid base derivative for a given target derivative
    
    # Arguments
    - `target_deriv::Dict{Symbol, Int}`: The target derivative dictionary
    - `computed::Set{String}`: Set of already computed derivative keys
    
    # Returns
    - `Tuple{Dict{Symbol, Int}, Symbol}`: A tuple containing the base derivative dictionary and the differentiating variable symbol
    
    # Notes
    - The function iterates through the variables in the target derivative and checks if the corresponding base derivative has been computed
    - If a valid base derivative is found, it returns the base derivative and the differentiating variable symbol
    - If no valid base derivative is found, returns (nothing, nothing)
=#
function find_valid_base_derivative(target_deriv::Dict{Symbol, Int}, computed::Set{String})::Tuple{Dict{Symbol, Int}, Symbol}
    for (var, order) in target_deriv
        if order > 0
            base_deriv = copy(target_deriv)
            base_deriv[var] = order - 1
            base_key = combine_derivative_strings(base_deriv)
            if base_key in computed
                return base_deriv, var
            end
        end
    end
    return nothing, nothing
end

#=  
    compute_third_derivatives_2D!(poly::Polynomial)
    Function to compute derivatives following two rules:
    Rule 1 : If the order of x is a multiple of 3 (including 0), split into two branches: one differentiating x and the other y.
    Rule 2 : If derivative order for y exists, then for all subsequent derivatives calculate w.r.t. y.
    
    # Arguments
    - `poly::Polynomial`: The polynomial structure containing the SLP and derivative information
    
    # Notes
    - If the polynomial does not have an SLP, an error is thrown
    - If the original function SLP range is missing, an error is thrown
    - The function computes derivatives for x and y according to the specified rules
=# 

function compute_third_derivatives_2D!(poly::Polynomial)
    empty!(INSTRUCTION_HASH)
    if poly.slp === nothing
        error("Polynomial must have an SLP to compute derivatives")
    end

    if !haskey(poly.slp.slp_ranges, "")
        error("Original function SLP range not found. The base SLP must be computed first.")
    end

    variables = get_variables(poly)


    if(length(variables) != 2)
        error("This partial derivative function is only defined for 2D polynomial")
    end
    
    to_compute = Vector{Dict{Symbol, Int}}()
    computed = Set{String}()

    push!(computed, "")

    for var in variables
        if get_max_order(poly, var) >= 1
            first_deriv = Dict{Symbol, Int}(var => 1)
            push!(to_compute, first_deriv)
        end
    end

    while !isempty(to_compute)
        current_deriv = popfirst!(to_compute)
        current_key = combine_derivative_strings(current_deriv)

        if current_key in computed
            continue
        end

        base_deriv, diff_var = find_valid_base_derivative(current_deriv, computed)

        if base_deriv !== nothing
            try 
                compute_derivative_slp!(poly, diff_var, base_deriv)
                push!(computed, current_key)

                x_order = get(current_deriv, :x, 0)
                y_order = get(current_deriv, :y, 0)
                
                if y_order > 0
                    new_deriv = copy(current_deriv)
                    current_y_order = get(new_deriv, :y, 0)
                    max_y_order = get_max_order(poly, :y)
                    
                    if current_y_order < max_y_order
                        new_deriv[:y] = current_y_order + 1
                        new_key = combine_derivative_strings(new_deriv)
                        
                        if !(new_key in computed) && !(new_deriv in to_compute)
                            push!(to_compute, new_deriv)
                        end
                    end
                elseif x_order % 3 == 0
                    for var in variables
                        new_deriv = copy(current_deriv)
                        current_order = get(new_deriv, var, 0)
                        max_order = get_max_order(poly, var)
                        
                        if current_order < max_order
                            new_deriv[var] = current_order + 1
                            new_key = combine_derivative_strings(new_deriv)

                            
                            if !(new_key in computed) && !(new_deriv in to_compute)
                                push!(to_compute, new_deriv)
                            end
                        end
                    end
                else
                    new_deriv = copy(current_deriv)
                    current_x_order = get(new_deriv, :x, 0)
                    max_x_order = get_max_order(poly, :x)
                    
                    if current_x_order < max_x_order
                        new_deriv[:x] = current_x_order + 1
                        new_key = combine_derivative_strings(new_deriv)
                        
                        if !(new_key in computed) && !(new_deriv in to_compute)
                            push!(to_compute, new_deriv)
                        end
                    end
                end
            catch e
                println("Error computing derivative for $current_key: $e")
            end
        else 
            push!(to_compute, current_deriv)
        end
    end
end




