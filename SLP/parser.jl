#= 
    @author: Thomas Chen
    @advisors: Yap Chee, Kai Hormann, Bingwei Zhang
    @date: 2025-06-17
    @description: parser
        The 2 main functions in this file are 
		   -- parse_expr()  which convert a string representation of a
		   		polynomial into an AST
           -- ast_to_slp() which convert an AST into an SLP
           
        The AST is a tree whose leaves are VarNodes and FloatNodes.
        The operation nodes are AddNode, MinusNode, MulNode, and PowNode.
    
		The SLP's have primitive instructions of the forms
			-- X+Y		(X,Y has the form of a, or a*(idx)
							where (idx) may be a variable or intermediate
							expression.
						(??? Do you allow negative constants ???)
			-- X-Y	
			-- X*Y
			-- X^n
						Constants are represented as a string
						that begins with "@"  (e.g., "@2.3" for 2.3)

		REMARK: We ought to give a BNF (Backus Normal Form)
			description of
				(1) Polynomial Expressions
				(2) SLP codes
			If we have these, then we can call tools like
				YACC or LEX to automatically generate these codes.
=# 

include("ASTTypes.jl")
include("myInterval.jl")
include("SLP.jl")

# Shortcut for zero interval
zero_interval = myInterval(0.0,0.0)


#=
    tokenize(input::String)::Vector{String}

    Parses the input String into tokens for AST Tree.

    # Arguments
    - `input::String`: The input string to be tokenized

    # Returns
    - `::Vector{String}`: A vector of tokens representing the input string
=#
function tokenize(input::String)::Vector{String}
    tokens = String[]
    i = 1
    len = length(input)

    while i <= len
        c = input[i]
        
        #Ignore white spaces
        if c in (' ', '\t', '\n')
            i += 1
            continue
        #Single character tokens
        elseif c in ('+', '-', '(', ')', '^', '*')
            push!(tokens, string(c))
            i += 1
        #Parse integer or float constants
        # Supports starting with ., ex. .7 = 0.7
        elseif isdigit(c) || c == '.' 
            start = i
            zero_prefix = (c == '.')
            has_decimal = zero_prefix
            i += 1
            if zero_prefix && (i > len || !isdigit(input[i]))
                error("Non-digit following decimal point")
            end

            while i <= len 
                c2 = input[i]
                if isdigit(c2)
                    i += 1
                elseif c2 == '.' && !has_decimal
                    has_decimal = true
                    i += 1

                    if i <= len && (!isdigit(input[i]) && !(input[i] in (' ', '\t', '\n'))) 
                        error("Non-digit following decimal point")
                    end
                elseif c2 == '.' && has_decimal
                    error("Only one decimal point allowed")
                else
                    break
                end
            end
            
            num = zero_prefix ? ("0" * input[start:i-1]) : input[start:i-1]
            num = endswith(num, ".") ? num[1:end-1] : num
            push!(tokens, num)

        #Parse Variable Name
         elseif isletter(c)
            token = string(c)
            i += 1
            while i <= len && isdigit(input[i])
                token *= string(input[i])
                i += 1
            end
            push!(tokens, token)
        else
            error("Unexpected character: $c")
        end 
    end 

    return tokens
end

#=
    Struct that records the state of the parser when parsing to AST

    tokens: List of tokens
    pos: current position(token index) the parser is working on
=#
mutable struct ParserState
    tokens::Vector{String}
    pos::Int
end

#=
    current_token(state::ParserState)::Union{String, Nothing}

    Retrieves the current_token that the parser is parsing

    # Arguments
    - `state::ParserState`: The current state of the parser

    # Returns
    - `::Union{String, Nothing}`: The current token being parsed. Returns `nothing` if at the end of the input.
=#
function current_token(state::ParserState)::Union{String, Nothing}
    if state.pos > length(state.tokens)
        return nothing
    end
    return state.tokens[state.pos]
end

#=
    consume_token!(state::ParserState)

    Moves the parser to the next token

    # Arguments
    - `state::ParserState`: The current state of the parser
=#
function consume_token!(state::ParserState)
    state.pos += 1
end


#=
    parse_expr(state::ParserState)::Expr

    The top level parsing function that parses an expr and separates expr into terms separated by + or - signs

    # Arguments
    - `state::ParserState`: The current state of the parser

    # Returns
    - `::Expr`: The parsed expression
=#
function parse_expr(state::ParserState)::Expr
    node  = parse_term(state)

    while true
        tok = current_token(state)
        if tok == "+"
            consume_token!(state)
            right = parse_term(state)
            node = AddNode(node, right)
        elseif tok == "-"
            consume_token!(state)
            right = parse_term(state)
            node = MinusNode(node, right)
        else 
            return node
        end
    end
end

#=
    parse_term(state::ParserState)::Expr

    The second level of parsing that takes a term and separates it into multiplication factors

    # Arguments
    - `state::ParserState`: The current state of the parser

    # Returns
    - `::Expr`: The parsed term expression
=#
function parse_term(state::ParserState)::Expr
    node = parse_factor(state)

    while true 
        tok = current_token(state)
        if tok == "*"
            consume_token!(state)
            right = parse_factor(state)
            node = MulNode(node, right)
        elseif tok !== nothing
            if tok == "(" || isletter(tok[1]) || isdigit(tok[1])
                right = parse_factor(state)
                node = MulNode(node, right)
            else
                return node
            end
        else
            return node
        end
    end
end

#=
    parse_factor(state::ParserState)::Expr

    The third level of parsing that takes a multiplication factor and retrieves its power

    # Arguments
    - `state::ParserState`: The current state of the parser

    # Returns
    - `::Expr`: The parsed factor expression
=#
function parse_factor(state::ParserState)::Expr
    node = parse_base(state)

    while true
        tok = current_token(state)
        if tok == "^"
            consume_token!(state)
            exponent_expr = parse_factor(state)
            if exponent_expr isa FloatNode
                node = PowNode(node, exponent_expr.value) 
            else
                error("Exponent must be a numeric literal, got $exponent_expr")
            end
        else 
            return node
        end
    end
end


#=
    parse_base(state::ParserState)::Expr

    The base level of parsing that takes a token and checks for its basic type(Variable, Float, or parenthesis) and converts them to ASTNode accordingly

    # Arguments
    - `state::ParserState`: The current state of the parser

    # Returns
    - `::Expr`: The parsed base expression
=#
function parse_base(state::ParserState)::Expr
    tok = current_token(state)
    if tok === nothing
        error("Unexpected end of input in parse_base")
    end


    if tok == "-"
        consume_token!(state)
        right = parse_base(state)
        return MinusNode(FloatNode(0.0), right)
    end
    
    if tok == "e"
        consume_token!(state)
        return FloatNode(exp(1.0))
    # variable
    elseif isletter(tok[1])
        consume_token!(state)
        return VarNode(Symbol(tok))
    # float64
    elseif isdigit(tok[1])
        consume_token!(state)
        return FloatNode(parse(Float64, tok))
    elseif tok == "("
        consume_token!(state)
        node = parse_expr(state)
        if current_token(state) != ")"
            error("Missing closing parenthesis")
        end
        consume_token!(state)
        return node
    else 
        error("Unexpected token in parse_base: $tok")
    end    
end

#=
    ast_to_slp(root::Expr)::SLP

    The function that takes an ast and converts it to a SLP. Mainly achieved using tree traversal technique. 

    # Arguments
    - `root::Expr`: The root of the AST tree

    # Returns
    - `::SLP`: The SLP representation of the AST
=#
function ast_to_slp(root::Expr)::SLP
    global INSTRUCTION_HASH
    # Clear the hash for fresh parsing
    empty!(INSTRUCTION_HASH)
    
    var_dict= Dict{Symbol,Int}()
    vars = Vector{Tuple{Int, Symbol}}()
    codeList = Vector{Tuple{Int, Symbol, operandType, operandType}}()

    nextVarIndex = Ref(-1)  
    nextTempIndex  = Ref(1)

    function get_index_for_var(var::Symbol)::Int
        if haskey(var_dict, var)
            return var_dict[var]
        else
            idx = nextVarIndex[]
            nextVarIndex[] -= 1
            var_dict[var] = idx
            push!(vars, (idx, var))
            return idx
        end
    end

    # Helper function to add instruction with hashset optimization
    function add_instruction_with_hash(op::Symbol, left_idx::Union{Int, String}, right_idx::Union{Int, String})::Int
        instruction_key = (op, left_idx, right_idx)
        
        # Check if instruction already exists
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
        
        # Create new instruction
        out_idx = nextTempIndex[]
        nextTempIndex[] += 1
        push!(codeList, (out_idx, op, left_idx, right_idx))
        
        # Add to hash
        INSTRUCTION_HASH[instruction_key] = out_idx
        if op == :+ || op == :*
            INSTRUCTION_HASH[(op, right_idx, left_idx)] = out_idx
        end
        
        return out_idx
    end

    function build(node::Expr)::Union{Int, String}
        if node isa VarNode
            return get_index_for_var(node.name)

        elseif node isa FloatNode
            return "@$(node.value)"

        elseif node isa AddNode
            left_idx = build(node.left)
            right_idx = build(node.right)
            
            # Constant folding optimizations
            if isa(left_idx, String) && startswith(left_idx, "@") && !contains(left_idx, "*")
                left_val = parse(Float64, left_idx[2:end])
                if left_val == 0.0
                    return right_idx  # 0 + x = x
                end
            end
        
            if isa(right_idx, String) && startswith(right_idx, "@") && !contains(right_idx, "*")
                right_val = parse(Float64, right_idx[2:end])
                if right_val == 0.0
                    return left_idx  
                end
            end
        
            if isa(left_idx, String) && isa(right_idx, String) && 
                startswith(left_idx, "@") && startswith(right_idx, "@") &&
                !contains(left_idx, "*") && !contains(right_idx, "*")
                left_val = parse(Float64, left_idx[2:end])
                right_val = parse(Float64, right_idx[2:end])
                return "@$(left_val + right_val)"
            end

            # Use hashset-optimized instruction addition
            return add_instruction_with_hash(:+, left_idx, right_idx)
    
        elseif node isa MinusNode
            left_idx = build(node.left)
            right_idx = build(node.right)

            if left_idx == right_idx 
                return "@0.0"
            end
            
            if isa(right_idx, String) && startswith(right_idx, "@") && !contains(right_idx, "*")
                right_val = parse(Float64, right_idx[2:end])
                if right_val == 0.0
                    return left_idx  # x - 0 = x
                end
            end

            if isa(left_idx, String) && isa(right_idx, String) && 
               startswith(left_idx, "@") && startswith(right_idx, "@") &&
               !contains(left_idx, "*") && !contains(right_idx, "*")
                left_val = parse(Float64, left_idx[2:end])
                right_val = parse(Float64, right_idx[2:end])
                return "@$(left_val - right_val)"
            end
            
            # Use hashset-optimized instruction addition
            return add_instruction_with_hash(:-, left_idx, right_idx)
    
        elseif node isa MulNode
            left_idx = build(node.left)
            right_idx = build(node.right)
            
            # Constant folding optimizations for multiplication
            if isa(left_idx, String) && startswith(left_idx, "@") && !contains(left_idx, "*")
                left_val = parse(Float64, left_idx[2:end])
                if left_val == 0.0
                    return "@0.0"  # 0 * x = 0
                elseif left_val == 1.0
                    return right_idx  # 1 * x = x
                end
            end
            
            if isa(right_idx, String) && startswith(right_idx, "@") && !contains(right_idx, "*")
                right_val = parse(Float64, right_idx[2:end])
                if right_val == 0.0
                    return "@0.0"  # x * 0 = 0
                elseif right_val == 1.0
                    return left_idx  # x * 1 = x
                end
            end

            if isa(left_idx, String) && isa(right_idx, String) && 
               startswith(left_idx, "@") && startswith(right_idx, "@") &&
               !contains(left_idx, "*") && !contains(right_idx, "*")
                left_val = parse(Float64, left_idx[2:end])
                right_val = parse(Float64, right_idx[2:end])
                return "@$(left_val * right_val)"
            end

            # Handle implicit multiplication cases
            if isa(left_idx, String) && startswith(left_idx, "@") && isa(right_idx, Int) && !contains(left_idx, "*")
                coeff = parse(Float64, left_idx[2:end])
                if coeff == 1.0
                    return right_idx 
                elseif coeff == 0.0
                    return "@0.0"     
                else
                    return "@$(coeff)*$(right_idx)"
                end
            end
            
            if isa(left_idx, Int) && isa(right_idx, String) && startswith(right_idx, "@") && !contains(right_idx, "*")
                coeff = parse(Float64, right_idx[2:end])
                if coeff == 1.0
                    return left_idx   
                elseif coeff == 0.0
                    return "@0.0"    
                else
                    return "@$(coeff)*$(left_idx)"
                end
            end
            
            # Use hashset-optimized instruction addition
            return add_instruction_with_hash(:*, left_idx, right_idx)
    
        elseif node isa PowNode
            n = node.exponent
            if isinteger(n)
                n = Int(n)
                if n == 0
                    return "@1.0" 
                elseif n == 1
                    return build(node.base)
                else
                    base_operand = build(node.base)
                    exp_operand = "@$(Float64(n))" 
                    
                    # Constant folding for exponentiation
                    if isa(base_operand, String) && startswith(base_operand, "@") && !contains(base_operand, "*")
                        base_val = parse(Float64, base_operand[2:end])
                        if base_val == 0.0
                            return "@0.0"  # 0^n = 0 (for n > 0)
                        elseif base_val == 1.0
                            return "@1.0"  # 1^n = 1
                        else
                            # Constant folding for power
                            return "@$(base_val^n)"
                        end
                    end
                    
                    # Use hashset-optimized instruction addition
                    return add_instruction_with_hash(:^, base_operand, exp_operand)
                end
            else
                error("Exponent is not an integer constant: $n")
            end
        else
            error("Unsupported AST node type: $(typeof(node))")
        end
    end

    build_result = build(root)
    
    slp_ranges = Dict{String, Tuple{Int, Int}}()
    
    if isempty(codeList)
        if !isempty(vars) && build_result isa Int
            push!(codeList, (1, :+, build_result, "@0.0"))
            slp_ranges[""] = (1, 1)
        elseif build_result isa String && startswith(build_result, "@")
            push!(codeList, (1, :+, "@0.0", build_result))
            slp_ranges[""] = (1, 1)
        else
            error("Unexpected empty codelist case")
        end
    else
        slp_ranges[""] = (1, length(codeList))
    end

    return SLP(vars, codeList, slp_ranges, Dict{Tuple{Int, Symbol}, Union{Int, String}}())
end

#=
    (Currently not used)
    parse_poly(poly::String)::SLP

    Top-level function that combines functions defined in this file and converts an input polynomial in string form to SLP form

    # Arguments
    - `poly::String`: The input polynomial in string form

    # Returns
    - `::SLP`: The SLP representation of the input polynomial
=#
function parse_poly(poly::String)::SLP
    if isempty(strip(poly))
        vars = Vector{Tuple{Int, Symbol}}()
        codeList = Vector{Tuple{Int, Symbol, operandType, operandType}}()
        slp_ranges = Dict{String, Tuple{Int, Int}}("" => (1, 0))
        return SLP(vars, codeList, slp_ranges, Dict{Tuple{Int, Symbol}, Union{Int, String}}())
    end

    tokenized_poly = tokenize(poly)
    state = ParserState(tokenized_poly, 1)
    poly_ast = parse_expr(state)
    slp = ast_to_slp(poly_ast)
    return slp
end
