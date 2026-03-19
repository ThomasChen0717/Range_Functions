#!/usr/bin/env julia

include("check_dep.jl")


using Printf
using BenchmarkTools
using Plots
using DataStructures
using Colors
using Dates
using NLsolve
using LinearAlgebra

include("SLP/myInterval.jl")
include("SLP/ASTTypes.jl")
include("SLP/SLP.jl")
include("SLP/utils.jl")
include("SLP/parser.jl")
include("SLP/polynomial.jl")
include("SLP/derivatives.jl")
include("SLP/eval.jl")


include("myBox.jl")
include("utils.jl")
include("tests.jl")

include("methods/lagrange3.jl")
include("methods/hermite4.jl")
include("methods/taylor4.jl")
include("methods/taylor3.jl")
include("methods/taylor2.jl")

#=
    Global configuration and shared state

    These variables cache derivative information and evaluation statistics that
    are reused by the different interval evaluation methods. They are populated
    and reset in helper routines (e.g. compute_all_derivatives!, reset_derivatives)
    and are treated as process-wide state during a run of the tool.

    - max_x, max_y          : maximum exponent of x and y in the current polynomial
    - total_degree          : total degree of the current polynomial
    - derivatives           : cached derivatives for generic methods
    - derivatives_taylor    : cached derivatives specialised for Taylor methods
    - derivatives_hermite   : cached derivatives specialised for Hermite methods
    - total_eval            : counter of total range evaluations performed
    - total_points          : counter of total sample points considered
    - tol                   : numerical tolerance used when filtering eigenvalues
    - factorial_cache       : cache of factorial values used by Taylor expansions
    - INSTRUCTION_HASH      : map used to deduplicate SLP instructions
=#
max_x, max_y = 0, 0
total_degree = 0
derivatives = Dict{Tuple{Float64, Float64}, Vector{Vector{Float64}}}()
derivatives_taylor = Dict{Tuple{Float64, Float64}, Vector{Vector{Float64}}}()
derivatives_hermite = Dict{Tuple{Float64, Float64}, Vector{Vector{Float64}}}()
total_eval = 0
total_points = 0
tol = 1e-10
factorial_cache = Vector{Int}()
INSTRUCTION_HASH = Dict{Tuple{Symbol, Union{Int, String}, Union{Int, String}}, Int}()

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 50
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10

#=
    main()

    Command-line entry point for the range functions analysis tool.

    Behaviour:
    - Parses ARGS to determine the command:
        * "intervals" : evaluate interval widths at several radii
        * "analyze"   : generate logD vs logr curves for two methods
        * "compare"   : generate a heatmap comparing two methods
    - Determines the list of input files (defaults to "input.txt" when omitted).
    - For each non-comment, non-empty line in each input file:
        * Uses parse_line to extract the polynomial string and x, y intervals.
        * Builds the corresponding myBox and its midpoint.
        * Dispatches to test_intervals_at_r, analyze, or compare_methods/
          create_visualization depending on the chosen command.
=#
function main()
    global max_x, max_y, total_degree

    cmd = ""
    input_files = String[]
    method1 = ""
    method2 = ""

    if isempty(ARGS)
        cmd = "intervals"
        input_files = ["input.txt"]
        println("No command specified. Using command 'intervals' and input file 'input.txt'.")
    elseif ARGS[1] == "compare"
        cmd = "compare"
        if length(ARGS) < 4
            error("Usage: compare METHOD1 METHOD2 D [input files...]")
        end
        method1 = ARGS[2]
        method2 = ARGS[3]
        d = parse(Int, ARGS[4])

        if length(ARGS) >= 5
            input_files = ARGS[5:end]
        else
            input_files = ["input.txt"]
        end
    elseif ARGS[1] == "analyze"
        cmd = "analyze"
        if length(ARGS) < 3
            error("Usage: analyze METHOD1 METHOD2 [input files...]")
        end
        method1 = ARGS[2]
        method2 = ARGS[3]
        if length(ARGS) >= 4
            input_files = ARGS[4:end]
        else
            input_files = ["input.txt"]
        end
    elseif ARGS[1] == "intervals"
        cmd = "intervals"
        if length(ARGS) >= 2
            input_files = ARGS[2:end]
        else
            input_files = ["input.txt"]
        end
    else
        cmd = "intervals"
        input_files = ARGS
        println("No command flag given. Interpreting arguments as input files.")
        println("Using default command 'intervals'.")
    end

    for input_file in input_files
        println("Reading input file: ", input_file)
        line_num = 0

        for line in eachline(input_file)
            line_num += 1
            s = strip(line)
            if isempty(s) || startswith(s, "#")
                continue
            end

            poly_str, varsDict = parse_line(String(s))

            if !haskey(varsDict, :x) || !haskey(varsDict, :y)
                println("Line ", line_num, ": missing x or y interval, skipping.")
                continue
            end

            x_val = varsDict[:x]
            y_val = varsDict[:y]

            if !(x_val isa myInterval && y_val isa myInterval)
                println("Line ", line_num, ": x and y must be intervals, skipping.")
                continue
            end

            B = myBox(x_val, y_val)
            midpoint_x = (x_val.lower + x_val.upper) / 2
            midpoint_y = (y_val.lower + y_val.upper) / 2

            println("Processing line ", line_num)
            println("Polynomial: ", poly_str)
            println("Box: x=", x_val, ", y=", y_val)

            if cmd == "intervals"
                test_intervals_at_r(poly_str, midpoint_x, midpoint_y)
            elseif cmd == "analyze"
                analyze(poly_str, midpoint_x, midpoint_y, method1, method2)
            elseif cmd == "compare"
                ratio_matrix = compare_methods(poly_str, B, method1, method2, d)
                p = create_visualization(ratio_matrix, B, poly_str, method1, method2)
                mkpath("imgs")
                poly_tag = replace(poly_str[1:min(20, length(poly_str))], r"[^a-zA-Z0-9]" => "_")
                filename = "imgs/compare_heatmap_$(poly_tag)_$(method1)_vs_$(method2)_x[$(B.x.lower)_$(B.x.upper)]_y[$(B.y.lower)_$(B.y.upper)].png"
                savefig(p, filename)
                println("Comparison plot saved to: ", filename)
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main() 
end
