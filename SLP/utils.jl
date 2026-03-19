
#= 
    @author: Thomas Chen
    @advisors: Yap Chee, Kai Hormann, Bingwei Zhang
    @date: 2025-06-17
    @description: Useful utility functions
=# 
using Printf

include("myInterval.jl")

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
    parse_line(line::String)::Tuple{String, Dict{Symbol, Union{Float64, myInterval}}, String}
     
    Function to parse the input file line by line and extract the desired parts.
    Supports variable substitutions that are either numbers (integers or floats)
    or intervals in the form "[a,b]".
    
    # Arguments
    - `line::String`: Input line from the file
    
    # Returns
    - `Tuple{String, Dict{Symbol, Union{Float64, myInterval}}, String}`: Tuple of polynomial, variables and description
=#
function parse_line(line::String)::Tuple{String, Dict{Symbol, Union{Float64, myInterval}}, String}
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

    description = ""
    if i <= length(tokens)
        description = join(tokens[i:end], ",")
        description = replace(description, "\"" => "")
    end

    return (poly, varsDict, description)
end

#=
    parse_derivative(key::String)::Dict{Symbol, Int}
    
    Parses a derivative string like "x^2y^1" into a dictionary {x: 2, y: 1}
    
    # Arguments
    - `key::String`: Input derivative string
    
    # Returns
    - `Dict{Symbol, Int}`: Dictionary mapping variable names to their derivative orders
=#
function parse_derivative(key::String)::Dict{Symbol, Int}
    if key == ""
        return Dict{Symbol, Int}()
    end
    
    result = Dict{Symbol, Int}()
    i = 1
    while i <= length(key)
        if isletter(key[i])
            var_name = string(key[i])
            i += 1
            
            # Check for ^power notation
            if i <= length(key) && key[i] == '^'
                i += 1  # skip ^
                power_str = ""
                while i <= length(key) && isdigit(key[i])
                    power_str *= key[i]
                    i += 1
                end
                power = parse(Int, power_str)
            else
                power = 1
            end
            
            result[Symbol(var_name)] = power
        else
            i += 1
        end
    end
    
    return result
end

#=
    derivative_string(var::Symbol, order::Int)::String
    
    Creates a derivative string for a single variable
    
    # Arguments
    - `var::Symbol`: Variable name
    - `order::Int`: Derivative order
    
    # Returns
    - `String`: Derivative string in the form "x^2"
=#
function derivative_string(var::Symbol, order::Int)::String
    if order == 0
        return ""
    elseif order == 1
        return string(var)
    else
        return string(var) * "^" * string(order)
    end
end

#=
    combine_derivative_strings(derivs::Dict{Symbol, Int})::String
    
    Combines multiple variable derivatives into a single string
    
    # Arguments
    - `derivs::Dict{Symbol, Int}`: Dictionary mapping variable names to their derivative orders
    
    # Returns
    - `String`: Combined derivative string in the form "x^2y^3"
=#
function combine_derivative_strings(derivs::Dict{Symbol, Int})::String
    if isempty(derivs)
        return ""
    end
    
    parts = String[]
    for (var, order) in sort(collect(derivs), by=x->string(x[1]))
        if order > 0
            push!(parts, derivative_string(var, order))
        end
    end
    
    return join(parts, "")
end

#=
    format_number(value::Float64, digits::Int, threshold::Float64=1e4)::String
     
    Format a Float64 number with specified precision digits.
    Uses scientific notation if the absolute value is larger than threshold.
    
    # Arguments
    - `value::Float64`: The number to format
    - `digits::Int`: Number of precision digits to display
    - `threshold::Float64`: Threshold above which to use scientific notation (default: 1e6)
    
    # Returns
    - `String`: Formatted number string
=#
function format_number(value::Float64, digits::Int,threshold::Int=6)::String
    abs_value = abs(value)
    
    if abs_value >= 10.0^threshold || (abs_value > 0 && abs_value < 10.0^-threshold)
        fmt = Printf.Format("%.$(digits)e")
        return Printf.format(fmt, value)
    else
        fmt = Printf.Format("%.$(digits)f")
        return Printf.format(fmt, value)
    end
end

#=
    format_number(interval::myInterval, digits::Int, threshold::Float64=1e6, centered::Bool=false)::String
     
    Format a myInterval with specified precision digits.
    Uses scientific notation if either bound is larger than threshold.
    When centered=true, displays as "center ± radius" format.
    
    # Arguments
    - `interval::myInterval`: The interval to format
    - `digits::Int`: Number of precision digits to display
    - `threshold::Int`: Threshold of order of magnitude above which to use scientific notation (default: 6)
    - `centered::Bool`: Whether to use centered form (center ± radius) or standard [lower, upper] form
    
    # Returns
    - `String`: Formatted interval string
=#
function format_number(interval::myInterval, digits::Int, threshold::Int=6, centered::Bool=false)::String
    lower_str = format_number(interval.lower, digits, threshold)
    upper_str = format_number(interval.upper, digits, threshold)
    if centered
        center = (interval.lower + interval.upper) / 2.0
        radius = (interval.upper - interval.lower) / 2.0
        
        center_str = format_number(center, digits, threshold)
        radius_str = format_number(radius, digits, threshold)
        
        return "$center_str ± $radius_str = [$lower_str, $upper_str]"
    else
        return "[$lower_str, $upper_str]"
    end
end

#=
    format_result(value::Union{Float64, myInterval}, digits::Int, threshold::Float64=1e6, centered::Bool=false)::String
     
    Generic function to format any numerical result (Float64 or myInterval).
    
    # Arguments
    - `value::Union{Float64, myInterval}`: The value to format
    - `digits::Int`: Number of precision digits to display
    - `threshold::Int`: Threshold of order of magnitude above which to use scientific notation (default: 6)
    - `centered::Bool`: Whether to use centered form for intervals
    
    # Returns
    - `String`: Formatted value string
=#
function format_result(value::Union{Float64, myInterval}, digits::Int, threshold::Int=6, centered::Bool=false)::String
    if value isa myInterval
        return format_number(value, digits, threshold, centered)
    else
        return format_number(value, digits, threshold)
    end
end



