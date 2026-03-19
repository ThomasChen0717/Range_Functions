#= 
    @author: Thomas Chen
    @advisors: Yap Chee, Kai Hormann, Bingwei Zhang
    @date: 2025-06-17
    @description: Polynomial Struct
=# 

include("ASTTypes.jl")
include("parser.jl")
include("SLP.jl")

#=
    Polynomial struct represents a polynomial abstractly
    
    Fields:
    - variables: Set of variable names (symbols) present in the polynomial
    - string_repr: String representation of the polynomial
    - max_orders: Dictionary mapping each variable to its maximum order in the polynomial
    - slp: SLP representation of the polynomial
=#
struct Polynomial
    variables::Set{Symbol}
    string_repr::String
    max_orders::Dict{Symbol, Int}
    total_degree::Int
    slp::Union{SLP, Nothing}
end

#=
    Constructor for Polynomial struct
    
    # Arguments
    - `poly_string::String`: String representation of the polynomial
    
    # Returns
    - `Polynomial`: Constructed Polynomial struct
=#
function Polynomial(poly_string::String)::Polynomial
    if isempty(strip(poly_string))
        vars = Vector{Tuple{Int, Symbol}}()
        codeList = Vector{Tuple{Int, Symbol, operandType,
				operandType}}()
        slp_ranges = Dict{String, Tuple{Int, Int}}("" => (1, 0))
        slp = SLP(vars, codeList, slp_ranges, Dict{Tuple{Int, Symbol},
				Union{Int, String}}())
        return Polynomial(Set{Symbol}(), poly_string,
				Dict{Symbol, Int}(), 0, slp)
    end

    tokenized_poly = tokenize(poly_string)
    state = ParserState(tokenized_poly, 1)
    poly_ast = parse_expr(state)
    slp = ast_to_slp(poly_ast)

    variables = Set{Symbol}()
    max_order = Dict{Symbol, Int}()
    total_degree = Ref(0)

    extract_variables_and_orders!(poly_ast, variables, max_order, total_degree)
    
    return Polynomial(variables, poly_string, max_order, total_degree[], slp)
end

#=
	Helper function to extract variables, their maximum orders, and total
		degree from AST
    
    # Arguments
    - `node::Expr`: The AST node to process
    - `variables::Set{Symbol}`: Set to store unique variable symbols
	- `max_order::Dict{Symbol, Int}`: Dictionary to store maximum order
	  		for each variable
	- `total_degree::Ref{Int}`: Reference to store the maximum total
	  		degree found
=#
function extract_variables_and_orders!(node::Expr,
		variables::Set{Symbol}, max_order::Dict{Symbol, Int},
		total_degree::Ref{Int})
    if node isa VarNode
        var_symbol = Symbol(node.name)
        push!(variables, var_symbol)
        max_order[var_symbol] = max(get(max_order, var_symbol, 0), 1)
        total_degree[] = max(total_degree[], 1)
    elseif node isa PowNode
        if node.base isa VarNode && node.exponent isa Number
            var_symbol = Symbol(node.base.name)
            push!(variables, var_symbol)
            order = Int(node.exponent)
            max_order[var_symbol] = max(get(max_order, var_symbol, 0), order)
            total_degree[] = max(total_degree[], order)
        else
            base_orders = Dict{Symbol, Int}()
            base_total_degree = Ref(0)
            extract_variables_and_orders!(node.base, variables, base_orders, base_total_degree)

            exponent_val = Int(node.exponent)
            for (var, base_order) in base_orders
                new_order = base_order * exponent_val
                max_order[var] = max(get(max_order, var, 0), new_order)
            end
            total_degree[] = max(total_degree[], base_total_degree[] * exponent_val)
        end
    elseif node isa MulNode
        left_orders = Dict{Symbol, Int}()
        right_orders = Dict{Symbol, Int}()
        left_total_degree = Ref(0)
        right_total_degree = Ref(0)
        
        extract_variables_and_orders!(node.left, variables, left_orders, left_total_degree)
        extract_variables_and_orders!(node.right, variables, right_orders, right_total_degree)
        
        all_vars = union(keys(left_orders), keys(right_orders))
        
        for var in all_vars
            left_degree = get(left_orders, var, 0)
            right_degree = get(right_orders, var, 0)
            total_var_degree = left_degree + right_degree
            max_order[var] = max(get(max_order, var, 0), total_var_degree)
        end
        
        # Total degree for multiplication is the sum of total degrees
        total_degree[] = max(total_degree[], left_total_degree[] + right_total_degree[])
    elseif node isa AddNode || node isa MinusNode
        left_total_degree = Ref(0)
        right_total_degree = Ref(0)
        
        extract_variables_and_orders!(node.left, variables, max_order, left_total_degree)
        extract_variables_and_orders!(node.right, variables, max_order, right_total_degree)
        
        # Total degree for addition/subtraction is the maximum of the two sides
        total_degree[] = max(total_degree[], left_total_degree[], right_total_degree[])
    end
end

# Convenience methods for the Polynomial struct
function get_variables(poly::Polynomial)::Set{Symbol}
    return poly.variables
end

function get_max_order(poly::Polynomial, var::Symbol)::Int
    return get(poly.max_orders, var, 0)
end

function get_total_degree(poly::Polynomial)::Int
    return poly.total_degree
end

function is_multivariate(poly::Polynomial)::Bool
    return length(poly.variables) > 1
end
