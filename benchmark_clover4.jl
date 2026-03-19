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

# Global configuration and shared state from range_funcs.jl
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

function run_clover4_benchmark()
    global max_x, max_y, total_degree

    poly_str = "50*x^10+(249*y^2-57)*x^8+(498*y^4-227*y^2-1)*x^6+(498*y^6-341*y^4-3*y^2+16)*x^4+(249*y^8-227*y^6-3*y^4-102*y^2-1)*x^2+50*y^10-57*y^8-y^6+16*y^4-y^2-1"


    B = myBox(myInterval(-1.2, 1.2), myInterval(-1.2, 1.2))

    println("Running benchmark for Clover4 polynomial:")
    println("Polynomial: ", poly_str)

    table_rows = NamedTuple{(:curve, :method, :min_ms, :median_ms, :mean_ms, :memory_mb, :efficacy), Tuple{String, String, Float64, Float64, Float64, Float64, Float64}}[]
    
    polynomial = Polynomial(poly_str)

    max_x, max_y = get_max_order(polynomial, :x), get_max_order(polynomial, :y)
    total_degree = get_total_degree(polynomial)

    description = "Clover4"
    box_count = 4096
    q = uniform_split(B, box_count)

    # methods = ["Lagrange3", "Lagrange3_shared", "Hermite4", "Hermite4_shared", "Taylor2", "Taylor3", "Taylor4"]
    methods = ["Taylor2"]

    for method in methods
        polynomial = Polynomial(poly_str)
        if method == "Lagrange3" || method == "Lagrange3_shared"
            compute_third_derivatives_2D!(polynomial)
        else
            compute_all_derivatives!(polynomial)
        end
        
        share = (method == "Lagrange3_shared" || method == "Hermite4_shared")
        eff_method = method
        if method == "Lagrange3_shared"
            eff_method = "Lagrange3"
        end 
        if method == "Hermite4_shared"
            eff_method = "Hermite4"
        end
        total_width = evaluate_boxes(q, polynomial, eff_method; sharing=share)
        efficacy = total_width / box_count

        benchmark_result = @benchmark evaluate_boxes($q, $polynomial, $eff_method; sharing=$share) teardown=(reset_derivatives())    

        min_ms    = minimum(benchmark_result.times) / 1e6
        median_ms = median(benchmark_result.times) / 1e6
        mean_ms   = mean(benchmark_result.times) / 1e6
        memory_mb = benchmark_result.memory / 1024^2

        println("\nResults for method: ", method)
        display(benchmark_result)
        println("Total width: $total_width")
        println("Efficacy: $efficacy")

        push!(table_rows, (; curve=description, method, min_ms, median_ms, mean_ms, memory_mb, efficacy))
    end

    println()
    println(rpad("Curve", 12), rpad("Method", 10), rpad("Min (ms)", 12),
            rpad("Median (ms)", 14), rpad("Mean (ms)", 12), rpad("Memory (MB)", 13), rpad("Efficacy", 12))

    for row in table_rows
        println(rpad(row.curve, 12),
                rpad(row.method, 10),
                @sprintf("%-11.4f", row.min_ms),
                @sprintf("%-13.4f", row.median_ms),
                @sprintf("%-11.4f", row.mean_ms),
                @sprintf("%-12.4f", row.memory_mb),
                @sprintf("%-12.4f", row.efficacy))
    end
    println()
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_clover4_benchmark()
end