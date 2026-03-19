#=
    monomial_string(ordX::Integer, ordY::Integer)::String

    Generates a monomial string based on the given orders of x and y.

    # Arguments
    - `ordX::Integer`: The order of x
    - `ordY::Integer`: The order of y

    # Returns
    - `::String`: The generated monomial string
=#
function monomial_string(ordX::Integer, ordY::Integer)::String
    parts = String[]

    if ordX == 1
        push!(parts, "x")
    elseif ordX > 1
        push!(parts, "x^" * string(ordX))
    end

    if ordY == 1
        push!(parts, "y")
    elseif ordY > 1
        push!(parts, "y^" * string(ordY))
    end

    return isempty(parts) ? "" : join(parts, "")
end

#=
    eval_slp(slp::SLP, x::Float64, y::Float64, order::String)::Float64

    Evaluates the SLP polynomial at a given point with a specified order.

    # Arguments
    - `slp::SLP`: The SLP polynomial
    - `x::Float64`: The x-coordinate
    - `y::Float64`: The y-coordinate
    - `order::String`: The order of the derivative

    # Returns
    - `::Float64`: The evaluated value
=#
function eval_slp(slp::SLP, x::Float64, y::Float64, order::String)::Float64
    global derivatives_taylor 
    point_key = (x, y)

    if haskey(derivatives_taylor , point_key)
        x_order, y_order = parse_monomial_key(order)
        
        cached_derivatives = derivatives_taylor[point_key]
        if x_order + 1 <= length(cached_derivatives) && y_order + 1 <= length(cached_derivatives[x_order + 1])
            return cached_derivatives[x_order + 1][y_order + 1]
        end
    end


    vars = Dict{Symbol,Union{Float64,myInterval}}()
    vars[:x] = x
    vars[:y] = y
    return evaluate_slp_range(slp, order, vars)
end

function interior!(a, b, c00,c10,c01,c20,c11,c02,c30,c21,c12,c03, r)
    # p_x = (3c30) x^2 + (2c20 + 2c21 y) x + (c10 + c11 y + c12 y^2)
    # p_y = ( c21) x^2 + (c11 + 2c12 y) x + (c01 + 2c02 y + 3c03 y^2)

    A = [3c30]                 # constant poly in y
    B = [2c20, 2c21]           # 2c20 + 2c21 y
    C = [c10,  c11,  c12]      # c10 + c11 y + c12 y^2

    D = [c21]
    E = [c11, 2c12]
    F = [c01, 2c02, 3c03]

    Res = resultant_quad_quad(A,B,C,D,E,F)
    ys  = real_roots_poly(Res)

    for y in ys
        if abs(y) > r + 1e-12
            continue
        end

        # solve p_x(x,y)=0 for x (quadratic in x)
        qa = 3c30
        qb = 2c20 + 2c21*y
        qc = c10 + c11*y + c12*y^2

        xs = Float64[]
        if abs(qa) < 1e-14
            if abs(qb) > 1e-14
                push!(xs, -qc/qb)
            end
        else
            disc = qb*qb - 4*qa*qc
            if disc ≥ -1e-12
                disc = max(disc, 0.0)
                sdisc = sqrt(disc)
                push!(xs, (-qb - sdisc)/(2*qa))
                push!(xs, (-qb + sdisc)/(2*qa))
            end
        end

        for x in xs
            if abs(x) > r + 1e-12
                continue
            end

            # filter: also satisfy p_y ≈ 0
            py = c21*x^2 + (c11 + 2c12*y)*x + (c01 + 2c02*y + 3c03*y^2)
            if abs(py) > 1e-8
                continue
            end

            p = c00 + c10*x + c01*y +
                c20*x^2 + c11*x*y + c02*y^2 +
                c30*x^3 + c21*x^2*y + c12*x*y^2 + c03*y^3

            a = min(a, p)
            b = max(b, p)
        end
    end

    return a, b
end

#=
real_roots_poly(coeffs::Vector{Float64}; tol_imag=1e-10)

Real roots of a polynomial with coefficients given low->high, computed
via eigenvalues of the companion matrix.
=#
function real_roots_poly(coeffs::Vector{Float64}; tol_imag=1e-10)
    n = length(coeffs) - 1
    n <= 0 && return Float64[]

    # convert to high->low
    a = reverse(coeffs)
    abs(a[1]) < 1e-14 && return Float64[]
    a ./= a[1]  # make monic

    C = zeros(Float64, n, n)
    if n > 1
        C[2:end, 1:end-1] .= I(n-1)
    end
    C[:, end] .= -a[2:end]

    vals = eigvals(C)

    roots = Float64[]
    for v in vals
        if abs(imag(v)) ≤ tol_imag
            push!(roots, real(v))
        end
    end
    return roots
end

#=
    get_point(B::myBox, i::Int, j::Int)::Tuple{Float64, Float64}

    Function to retrieve a specific point from the box region.

    # Arguments
    - `B`: The box region
    - `i`: The index of the point in the x direction
    - `j`: The index of the point in the y direction

    # Returns
    - `Tuple{Float64, Float64}`: A tuple containing the coordinates of the point
=#
function get_point(B::myBox, i::Int, j::Int)::Tuple{Float64, Float64}
    return B.pts_matrix[i, j]
end

#= 
    (For testing)
    uniform_split(b::myBox, dim::Int)::PriorityQueue{myBox, Int}

    Splits a box into `dim` subboxes of equal size and enqueues them with
    depth-based priority.

    # Arguments
    - `b::myBox`: The box to split
    - `dim::Int`: Target number of subboxes (must be ≥ 1)

    # Returns
    - `PriorityQueue{myBox, Int}`: Priority queue containing the split boxes
=#
function uniform_split(b::myBox, dim::Int)::PriorityQueue{myBox, Int}
    box_queue = PriorityQueue{myBox, Int}()
    enqueue!(box_queue, b, b.depth)

    while length(box_queue) < dim
        box_to_split = dequeue!(box_queue)
        children = split_box(box_to_split)

        for child in children 
            enqueue!(box_queue, child, child.depth)
        end
    end
    return box_queue
end

#= 
    (For testing)
    evaluate_boxes(q::PriorityQueue{myBox, Int}, poly::Polynomial, method::String; sharing::Bool = true)

    Evaluates all boxes in the priority queue using the specified interval
    method (Lagrange3, Taylor2/3/4, Hermite4) and accumulates the total width
    of the resulting intervals.

    # Arguments
    - `q::PriorityQueue{myBox, Int}`: Priority queue containing the boxes
    - `poly::Polynomial`: Polynomial to be evaluated
    - `method::String`: Name of the interval method to use
    - `sharing::Bool`: Whether to share derivative caches between boxes

    # Returns
    - `Float64`: Sum of interval widths over all boxes
=#
function evaluate_boxes(q::PriorityQueue{myBox, Int}, poly::Polynomial, method::String; sharing::Bool = true)
    global total_degree

    total_width = 0.0
    for b in keys(q)
        if !sharing
            reset_derivatives()
        end
        if method == "Taylor4"
            range = taylor_interpolation4(poly, b, total_degree; sharing=sharing)
        elseif method == "Taylor3"
            range = taylor_interpolation3(poly, b, total_degree; sharing=sharing)
        elseif method == "Taylor2"
            range = taylor_interpolation2(poly, b, total_degree; sharing=sharing)
        elseif method == "Hermite4"
            range = hermite4(poly, b, total_degree; sharing=sharing)
        elseif method == "Lagrange3"
            range = Lagrange3(poly, b, total_degree; sharing=sharing)
        end
        curr_width = range.upper - range.lower
        total_width += curr_width
    end

    return total_width
end

function reset_derivatives()
    global derivatives, derivatives_taylor, derivatives_hermite
    derivatives = Dict{Tuple{Float64, Float64}, Vector{Vector{Float64}}}()
    derivatives_taylor = Dict{Tuple{Float64, Float64}, Vector{Vector{Float64}}}()
    derivatives_hermite = Dict{Tuple{Float64, Float64}, Vector{Vector{Float64}}}()
end


function printf(f)
    if abs(f) > 1e10
        @printf(stdout,"%.3e",f)
    else 
        @printf(stdout,"%.3f",f)
    end
end

function get_width(I::myInterval)
    return I.upper - I.lower
end

function get_width(I::Tuple{Float64, Float64})
    return I[2] - I[1]
end

#=
    writeTitle(title::String, file::IOStream)
     
    Function to write file name for each file read for visual clarity
    
    # Arguments
    - `title::String`: Title of the file to be written
    - `file::IOStream`: File stream to write the title
=#
function writeTitle(title::String, file::IOStream)
    n = length(title)+1
    write(file, string(repeat("-", n), "\n"))
    write(file, "$(title):\n")
    write(file, string(repeat("-", n), "\n"))
end

#=
    parse_line(line::String)::Tuple{String, Dict{Symbol, Union{Float64, myInterval}}}
     
    Function to parse the input file line by line and extract the desired parts.
    Supports variable substitutions that are either numbers (integers or floats)
    or intervals in the form "[a,b]".
    
    # Arguments
    - `line::String`: Input line from the file
    
    # Returns
    - `Tuple{String, Dict{Symbol, Union{Float64, myInterval}}}`: Tuple of polynomial, variables
=#
function parse_line(line::String)::Tuple{String, Dict{Symbol, Union{Float64, myInterval}}}
    tokens = String[]
    current = IOBuffer()
    in_brackets = false
    in_quotes = false

    for c in line
        if c == '"'
            in_quotes = !in_quotes
            write(current, c)
        elseif c == '[' && !in_quotes
            in_brackets = true
            write(current, c)
        elseif c == ']' && !in_quotes
            in_brackets = false
            write(current, c)
        elseif c == ',' && !in_brackets && !in_quotes
            token = String(take!(current))
            push!(tokens, strip(token))
        else
            write(current, c)
        end
    end

    token = String(take!(current))
    if !isempty(strip(token))
        push!(tokens, strip(token))
    end

    if isempty(tokens)
        error("No tokens found in input line")
    end

    poly = tokens[1]

    varsDict = Dict{Symbol,Union{Float64,myInterval}}()
    i = 2

    while i <= length(tokens) && !startswith(tokens[i], "[")
        if occursin("=", tokens[i])
            parts = split(tokens[i], "=")
            var = Symbol(strip(parts[1]))
            valstr = strip(parts[2])

            if startswith(valstr, "[") && endswith(valstr, "]")
                inner = valstr[2:end-1]
                splitVals = split(inner, ",")
                if length(splitVals) != 2
                    error("Interval for variable $(var) must have two endpoints")
                end
                a = parse(Float64, strip(splitVals[1]))
                b = parse(Float64, strip(splitVals[2]))
                varsDict[var] = myInterval(a, b)
            else
                numVal = parse(Float64, valstr)
                varsDict[var] = numVal
            end
        end
        i += 1
    end

    return (poly, varsDict)
end
