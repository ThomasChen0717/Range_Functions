#=
    critical_point_range(box::myBox, poly::Polynomial)::myInterval
    
    Compute the exact range of f over `box` by locating true extremal points.
    Method:
    1. Interior critical points: solve ∇f = 0
    2. Boundary extrema: on each edge solve df/ds = 0 (derivative along the edge direction)
    3. Corner points: evaluate the four corners
    
    Arguments:
    - `box::myBox`: computation region
    - `poly::Polynomial`: polynomial to evaluate
    
    Returns:
    - `myInterval`: exact range of f(B)
=#
function critical_point_range(box::myBox, poly::Polynomial)::myInterval
    x_1, x_3 = box.x.lower, box.x.upper
    y_1, y_3 = box.y.lower, box.y.upper
    
    candidate_values = Float64[]
    
    vars = Dict{Symbol, Union{Float64, myInterval}}()
    
    # ==================== 1. Evaluate the four corner points ====================
    corners = [(x_1, y_1), (x_1, y_3), (x_3, y_1), (x_3, y_3)]
    for (x, y) in corners
        vars[:x], vars[:y] = x, y
        val = evaluate_slp_range(poly.slp, "", vars)
        push!(candidate_values, val)
    end
    
    # ==================== 2. Interior critical points: ∇f = 0 ====================
    function gradient_func(vars_vec::Vector{Float64})
        x, y = vars_vec[1], vars_vec[2]
        
        var_dict_x = Dict{Symbol, Union{Float64, myInterval}}(:x => x, :y => y)
        var_dict_y = Dict{Symbol, Union{Float64, myInterval}}(:x => x, :y => y)
        
        fx = evaluate_slp_range(poly.slp, "x", var_dict_x)
        fy = evaluate_slp_range(poly.slp, "y", var_dict_y)
        
        return [fx, fy]
    end
    
    # Search interior critical points with multiple initial guesses
    n_search = 5
    for i in 1:n_search
        for j in 1:n_search
            x0 = x_1 + (x_3 - x_1) * (i / (n_search + 1))
            y0 = y_1 + (y_3 - y_1) * (j / (n_search + 1))
            
            try
                result = nlsolve(gradient_func, [x0, y0], method=:trust_region, ftol=1e-10)
                
                if result.x_converged || result.f_converged
                    x_sol, y_sol = result.zero[1], result.zero[2]
                    
                    # Strict check: must be inside the box (excluding boundary)
                    if x_1 < x_sol < x_3 && y_1 < y_sol < y_3
                        vars[:x], vars[:y] = x_sol, y_sol
                        val = evaluate_slp_range(poly.slp, "", vars)
                        push!(candidate_values, val)
                    end
                end
            catch
            end
        end
    end
    
    # ==================== 3. Bottom edge boundary extrema: y = y_1, solve df/dx = 0 ====================
    # On the bottom edge: extrema of f(x, y_1) along x satisfy df/dx = 0
    function bottom_edge_deriv(x_vec::Vector{Float64})
        x = x_vec[1]
        var_dict = Dict{Symbol, Union{Float64, myInterval}}(:x => x, :y => y_1)
        fx = evaluate_slp_range(poly.slp, "x", var_dict)
        return [fx]
    end
    
    try
        for x_init in range(x_1 + (x_3-x_1)*0.1, x_3 - (x_3-x_1)*0.1, length=3)
            result = nlsolve(bottom_edge_deriv, [x_init], method=:trust_region, ftol=1e-10)
            if result.x_converged || result.f_converged
                x_sol = result.zero[1]
                if x_1 < x_sol < x_3
                    vars[:x], vars[:y] = x_sol, y_1
                    val = evaluate_slp_range(poly.slp, "", vars)
                    push!(candidate_values, val)
                end
            end
        end
    catch
    end
    
    # ==================== 4. Top edge boundary extrema: y = y_3, solve df/dx = 0 ====================
    function top_edge_deriv(x_vec::Vector{Float64})
        x = x_vec[1]
        var_dict = Dict{Symbol, Union{Float64, myInterval}}(:x => x, :y => y_3)
        fx = evaluate_slp_range(poly.slp, "x", var_dict)
        return [fx]
    end
    
    try
        for x_init in range(x_1 + (x_3-x_1)*0.1, x_3 - (x_3-x_1)*0.1, length=3)
            result = nlsolve(top_edge_deriv, [x_init], method=:trust_region, ftol=1e-10)
            if result.x_converged || result.f_converged
                x_sol = result.zero[1]
                if x_1 < x_sol < x_3
                    vars[:x], vars[:y] = x_sol, y_3
                    val = evaluate_slp_range(poly.slp, "", vars)
                    push!(candidate_values, val)
                end
            end
        end
    catch
    end
    
    # ==================== 5. Left edge boundary extrema: x = x_1, solve df/dy = 0 ====================
    function left_edge_deriv(y_vec::Vector{Float64})
        y = y_vec[1]
        var_dict = Dict{Symbol, Union{Float64, myInterval}}(:x => x_1, :y => y)
        fy = evaluate_slp_range(poly.slp, "y", var_dict)
        return [fy]
    end
    
    try
        for y_init in range(y_1 + (y_3-y_1)*0.1, y_3 - (y_3-y_1)*0.1, length=3)
            result = nlsolve(left_edge_deriv, [y_init], method=:trust_region, ftol=1e-10)
            if result.x_converged || result.f_converged
                y_sol = result.zero[1]
                if y_1 < y_sol < y_3
                    vars[:x], vars[:y] = x_1, y_sol
                    val = evaluate_slp_range(poly.slp, "", vars)
                    push!(candidate_values, val)
                end
            end
        end
    catch
    end
    
    # ==================== 6. Right edge boundary extrema: x = x_3, solve df/dy = 0 ====================
    function right_edge_deriv(y_vec::Vector{Float64})
        y = y_vec[1]
        var_dict = Dict{Symbol, Union{Float64, myInterval}}(:x => x_3, :y => y)
        fy = evaluate_slp_range(poly.slp, "y", var_dict)
        return [fy]
    end
    
    try
        for y_init in range(y_1 + (y_3-y_1)*0.1, y_3 - (y_3-y_1)*0.1, length=3)
            result = nlsolve(right_edge_deriv, [y_init], method=:trust_region, ftol=1e-10)
            if result.x_converged || result.f_converged
                y_sol = result.zero[1]
                if y_1 < y_sol < y_3
                    vars[:x], vars[:y] = x_3, y_sol
                    val = evaluate_slp_range(poly.slp, "", vars)
                    push!(candidate_values, val)
                end
            end
        end
    catch
    end
    
    # ==================== Compute minimum and maximum ====================
    if isempty(candidate_values)
        # Fallback: if no critical points are found, return corner range
        return myInterval(
            minimum([evaluate_slp_range(poly.slp, "", Dict{Symbol, Union{Float64, myInterval}}(:x => x, :y => y)) for (x, y) in corners]),
            maximum([evaluate_slp_range(poly.slp, "", Dict{Symbol, Union{Float64, myInterval}}(:x => x, :y => y)) for (x, y) in corners])
        )
    end
    
    min_val = minimum(candidate_values)
    max_val = maximum(candidate_values)
    
    return myInterval(min_val, max_val)
end

#=
    uniform_sample_range(box::myBox, poly::Polynomial, num_samples::Int=100)::myInterval
    
    Compute the exact range f(B) using the critical point method rather than simple sampling.
    
    Arguments:
    - `box::myBox`
    - `poly::Polynomial`
    - `num_samples::Int`: deprecated, kept for compatibility
    
    Returns:
    - `myInterval`
=#
function uniform_sample_range(box::myBox, poly::Polynomial, num_samples::Int=100)::myInterval
    return critical_point_range(box, poly)
end

# Keep the old sampling method as a fallback
function uniform_sample_range_old(box::myBox, poly::Polynomial, num_samples::Int=100)::myInterval

    x_min = box.x.lower
    x_max = box.x.upper
    y_min = box.y.lower
    y_max = box.y.upper
    
    samples_per_dim = Int(ceil(sqrt(num_samples)))
    
    x_step = (x_max - x_min) / (samples_per_dim - 1)
    y_step = (y_max - y_min) / (samples_per_dim - 1)
    
    min_val = Inf
    max_val = -Inf
    
    vars = Dict{Symbol,Union{Float64,myInterval}}()


    for i in 0:(samples_per_dim-1)
        for j in 0:(samples_per_dim-1)
            x = x_min + i * x_step
            y = y_min + j * y_step

            vars[:x] = x
            vars[:y] = y
            val = evaluate_slp_range(poly.slp, "", vars)
        
            if val < min_val
                min_val = val
            end
            if val > max_val
                max_val = val
            end
        end
    end
    
    return myInterval(min_val, max_val)
end


# Compute the logarithmic distance between a method's interval and the
# reference (real) interval:
#     logD = log10(max(|method_lower - real_lower|,
#                      |method_upper - real_upper|)).
# If the intervals coincide exactly, returns -16.0 as a floor value.
function compute_logD(method_interval::myInterval, real_interval::myInterval)::Float64
    lower_diff = abs(method_interval.lower - real_interval.lower)
    upper_diff = abs(method_interval.upper - real_interval.upper)
    max_diff = max(lower_diff, upper_diff)
    if max_diff == 0.0
        return -16.0
    end
    
    return log10(max_diff)
end

# Create axis-aligned boxes of varying radius around a midpoint.
# Each radius r generates a box
#     [midpoint_x - r, midpoint_x + r] × [midpoint_y - r, midpoint_y + r].
function create_boxes_with_radius(midpoint_x::Float64, midpoint_y::Float64, radii::Vector{Float64})::Vector{myBox}
    boxes = myBox[]
    
    for r in radii
        box = myBox(
            myInterval(midpoint_x - r, midpoint_x + r),
            myInterval(midpoint_y - r, midpoint_y + r)
        )
        push!(boxes, box)
    end
    
    return boxes
end

# Compare two enclosure methods by plotting logD vs log r for shrinking
# square boxes around (midpoint_x, midpoint_y).  Uses method1 and method2
# as in compare_methods (Lagrange3, Taylor*, Hermite4).
function analyze(poly_str::String, midpoint_x::Float64, midpoint_y::Float64, method1::String, method2::String)
    poly_lagrange = Polynomial(poly_str)
    poly_taylor = Polynomial(poly_str)

    global max_x, max_y, total_degree
    max_x, max_y = get_max_order(poly_lagrange, :x), get_max_order(poly_lagrange, :y)
    total_degree = get_total_degree(poly_lagrange)

    if method1 == "Lagrange3" || method2 == "Lagrange3"
        compute_third_derivatives_2D!(poly_lagrange)
    end

    if method1 in ("Taylor4", "Taylor3", "Taylor2", "Hermite4") || method2 in ("Taylor4", "Taylor3", "Taylor2", "Hermite4")
        compute_all_derivatives!(poly_taylor)
    end

    radii = [10^x for x in range(log10(0.0000001), log10(1.0), length=1000)]
    log_radii = log10.(radii)

    boxes = create_boxes_with_radius(midpoint_x, midpoint_y, radii)

    method1_logD = Float64[]
    method2_logD = Float64[]

    function method_interval(poly_lagrange::Polynomial, poly_taylor::Polynomial, box::myBox, method::String)::myInterval
        if method == "Lagrange3"
            return Lagrange3(poly_lagrange, box, total_degree)
        elseif method == "Hermite4"
            return hermite4(poly_taylor, box, total_degree)
        elseif method == "Taylor4"
            return taylor_interpolation4(poly_taylor, box, total_degree)
        elseif method == "Taylor3"
            return taylor_interpolation3(poly_taylor, box, total_degree)
        elseif method == "Taylor2"
            return taylor_interpolation2(poly_taylor, box, total_degree)
        else
            error("Unknown method: $method")
        end
    end

    for (i, (box, r)) in enumerate(zip(boxes, radii))
        if i % 10 == 0
            println("Processing box $i/$(length(boxes)), radius = $r")
        end
        
        try
            real_interval = uniform_sample_range(box, poly_lagrange)

            reset_derivatives()
            fill!(box.FB, nothing)
            fill!(box.QB, [])
    
            interval1 = method_interval(poly_lagrange, poly_taylor, box, method1)

            reset_derivatives()
            fill!(box.FB, nothing)
            fill!(box.QB, [])
            
            interval2 = method_interval(poly_lagrange, poly_taylor, box, method2)

            reset_derivatives()
            fill!(box.FB, nothing)
            fill!(box.QB, [])
            
            logd1 = compute_logD(interval1, real_interval)
            logd2 = compute_logD(interval2, real_interval)
            
            push!(method1_logD, logd1)
            push!(method2_logD, logd2)
            
        catch e
            println("Error processing box with radius $r")
            push!(method1_logD, NaN)
            push!(method2_logD, NaN)
        end
    end
    
    valid_indices = .!isnan.(method1_logD) .& .!isnan.(method2_logD)
    if !any(valid_indices)
        println("Analysis produced no valid data points; skipping plot.")
        return nothing, log_radii, method1_logD, method2_logD
    end
    log_radii_clean = log_radii[valid_indices]
    method1_logD_clean = method1_logD[valid_indices]
    method2_logD_clean = method2_logD[valid_indices]
    
    p = plot(log_radii_clean, method1_logD_clean, 
             label=method1, 
             linewidth=0.5, 
             marker=:circle, 
             markersize=0.5,
             markerstrokewidth=0.5,
             color=:blue,
             markerstrokecolor=:darkblue,
             alpha=0.8,
             xlabel="log₁₀(r)", 
             ylabel="log₁₀(D)",
             title="logD vs logr Analysis\nPolynomial: $(poly_str[1:min(50, length(poly_str))])...",
             size=(800, 600),
             dpi=300,
             legend=:topright)
    
    plot!(p, log_radii_clean, method2_logD_clean, 
          label=method2, 
          linewidth=0.5, 
          marker=:square, 
          markersize=0.5,
          markerstrokewidth=0.5,
          color=:red,
          markerstrokecolor=:darkred,
          alpha=0.8,
          linestyle=:solid)

    x_range = extrema(log_radii_clean)
    y_intercept = minimum([minimum(method1_logD_clean), minimum(method2_logD_clean)]) - 1
    reference_y = [y_intercept + 3 * (x - x_range[1]) for x in log_radii_clean]
    
    plot!(p, log_radii_clean, reference_y,
          label="Reference (slope=3)",
          linewidth=2,
          linestyle=:dash,
          color=:green,
          alpha=0.7)
    
    plot!(p, grid=true, gridwidth=1, gridcolor=:gray, gridalpha=0.3)
    
    println("Analysis complete. Generated plot with $(length(log_radii_clean)) valid data points.")

    mkpath("imgs")

    filename = "imgs/logD_vs_logr_analysis_$(replace(poly_str[1:min(20, length(poly_str))], r"[^a-zA-Z0-9]" => "_"))_$(method1)_vs_$(method2)_x$(midpoint_x)_y$(midpoint_y).png"
    savefig(p, filename)
    
    println("Plot saved to: $filename")
    
    return p, log_radii_clean, method1_logD_clean, method2_logD_clean
end


# Evaluate several interval methods on a small list of radii and print
# their enclosures alongside a Monte‑Carlo reference interval.
function test_intervals_at_r(poly_str::String, midpoint_x::Float64, midpoint_y::Float64; radii::Vector{Float64}=[0.2, 0.1, 0.05])
    poly_lagrange = Polynomial(poly_str)
    poly_taylor = Polynomial(poly_str)
    poly_hermite = Polynomial(poly_str)

    global max_x, max_y, total_degree
    max_x, max_y = get_max_order(poly_lagrange, :x), get_max_order(poly_lagrange, :y)
    total_degree = get_total_degree(poly_lagrange)

    compute_third_derivatives_2D!(poly_lagrange)
    compute_all_derivatives!(poly_taylor)
    compute_all_derivatives!(poly_hermite)

    results = Tuple{Float64, myInterval, myInterval, myInterval}[]

    for r in radii
        box = myBox(
            myInterval(midpoint_x - r, midpoint_x + r),
            myInterval(midpoint_y - r, midpoint_y + r)
        )

        real_interval = uniform_sample_range(box, poly_lagrange)

        reset_derivatives()
        fill!(box.FB, nothing)
        fill!(box.QB, [])
        poly_lagrange = Polynomial(poly_str)
        compute_third_derivatives_2D!(poly_lagrange)

        # lagrange_interval = lagrange_range_function(poly_lagrange, box, total_degree, "S", 6)
        lagrange_interval = Lagrange3(poly_lagrange, box, total_degree)

        reset_derivatives()
        fill!(box.FB, nothing)
        fill!(box.QB, [])

        taylor_interval_3 = taylor_interpolation3(poly_taylor, box, total_degree)

        reset_derivatives()
        fill!(box.FB, nothing)
        fill!(box.QB, [])

        poly_taylor = Polynomial(poly_str)
        compute_all_derivatives!(poly_taylor)


        taylor_interval_2 = taylor_interpolation2(poly_taylor, box, total_degree)
       
        reset_derivatives()
        fill!(box.FB, nothing)
        fill!(box.QB, [])
        poly_taylor = Polynomial(poly_str)
        compute_all_derivatives!(poly_taylor)


        taylor_interval_4 = taylor_interpolation4(poly_taylor, box, total_degree)
       
        reset_derivatives()
        fill!(box.FB, nothing)
        fill!(box.QB, [])

        hermite_interval = hermite4(poly_hermite, box, total_degree)
        # hermite_interval = hermite4(poly_str, box_center, total_degree)
        println("r = $(r)")
        println("real_interval = $(real_interval)")
        println("lagrange_interval = $(lagrange_interval)")
        println("taylor_interval_2 = $(taylor_interval_2)")
        println("taylor_interval_3 = $(taylor_interval_3)")
        println("taylor_interval_4 = $(taylor_interval_4)")
        println("hermite_interval = $(hermite_interval)")

        push!(results, (r, real_interval, lagrange_interval, taylor_interval_3))
    end

    return results
end

function compare_methods(poly_str::String, box::myBox, method1::String, method2::String, dim::Int = 1024)
    #=
    Compare two interval methods on subdivided boxes.
    `dim` controls the total number of subboxes (e.g. 32x32 = 1024).
    Returns a matrix of ratios: width(method1) / width(method2).
    =#
    poly_lagrange = Polynomial(poly_str)
    poly_taylor = Polynomial(poly_str)
    
    global max_x, max_y, total_degree
    max_x, max_y = get_max_order(poly_lagrange, :x), get_max_order(poly_lagrange, :y)
    total_degree = get_total_degree(poly_lagrange)

    if method1 == "Lagrange3" || method2 == "Lagrange3"
        compute_third_derivatives_2D!(poly_lagrange)
    end

    if method1 in ("Taylor4", "Taylor3", "Taylor2", "Hermite4") || method2 in ("Taylor4", "Taylor3", "Taylor2", "Hermite4")
        compute_all_derivatives!(poly_taylor)
    end

    function method_width(poly_lagrange::Polynomial, poly_taylor::Polynomial, subbox::myBox, method::String)
        if method == "Lagrange3"
            interval = Lagrange3(poly_lagrange, subbox, total_degree)
        elseif method == "Hermite4"
            interval = hermite4(poly_taylor, subbox, total_degree)
        elseif method == "Taylor4"
            interval = taylor_interpolation4(poly_taylor, subbox, total_degree)
        elseif method == "Taylor3"
            interval = taylor_interpolation3(poly_taylor, subbox, total_degree)
        elseif method == "Taylor2"
            interval = taylor_interpolation2(poly_taylor, subbox, total_degree)
        else
            error("Unknown method: $method, must be one of Lagrange3, Hermite4, Taylor4, Taylor3, Taylor2")
        end
        return get_width(interval)
    end

    side = round(Int, sqrt(dim))
    ratio_matrix = zeros(Float64, side, side)
    
    x_step = (box.x.upper - box.x.lower) / side
    y_step = (box.y.upper - box.y.lower) / side
    
    println("Computing ranges for $(side)x$(side) = $(side*side) subboxes using methods $method1 vs $method2...")
    
    for i in 1:side
        for j in 1:side
            x_lower = box.x.lower + (j-1) * x_step
            x_upper = box.x.lower + j * x_step
            y_lower = box.y.lower + (i-1) * y_step
            y_upper = box.y.lower + i * y_step
            
            subbox = myBox(myInterval(x_lower, x_upper), myInterval(y_lower, y_upper))

            box_center = ((x_lower + x_upper) / 2, (y_lower + y_upper) / 2, (x_upper - x_lower) / 2)
            
            println("Processing box at ($i, $j) with center $box_center")

            if (i-1)*side + j % 100 == 0
                println("Processed $((i-1)*side + j)/$(side*side) boxes")
            end
            
            try
                width1 = method_width(poly_lagrange, poly_taylor, subbox, method1)
                reset_derivatives()

                width2 = method_width(poly_lagrange, poly_taylor, subbox, method2)
                reset_derivatives()
                
                if width2 > 0
                    ratio = width1 / width2
                else
                    ratio = 1.0
                end
                
                ratio_matrix[i, j] = ratio
                
            catch e
                println("Error processing box at ($i, $j) for methods $method1 vs $method2: $e")
                ratio_matrix[i, j] = 1.0 
            end
        end
    end
    
    return ratio_matrix
end
# Visualize width ratios as a heatmap over box, with green where method1
# is tighter (ratio ≤ 1) and red where method2 is tighter (ratio > 1),
# and overlay the zero level set of the polynomial.
function create_visualization(ratio_matrix::Matrix{Float64}, box::myBox, poly_str::String, method1::String, method2::String)
    color_matrix = zeros(size(ratio_matrix))
    
    for i in 1:size(ratio_matrix, 1)
        for j in 1:size(ratio_matrix, 2)
            ratio = ratio_matrix[i, j]
            
            if ratio > 1.0
                # Red side: Taylor is better (ratio > 1)
                color_matrix[i, j] = log10(ratio)
            else
                # Green side: Lagrange is better (ratio ≤ 1)
                color_matrix[i, j] = -abs(log10(ratio))
            end
        end
    end
    
    x_range = range(box.x.lower, box.x.upper, length=size(ratio_matrix, 2))
    y_range = range(box.y.lower, box.y.upper, length=size(ratio_matrix, 1))
    
    p = heatmap(x_range, y_range, color_matrix,
               color=cgrad([:darkgreen, :forestgreen, :limegreen, :yellow, :orange, :red, :darkred], [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]),
               aspect_ratio=:equal,
               title="$method1 vs $method2 Tightness Comparison\n$(poly_str[1:min(50, length(poly_str))])...",
               xlabel="",
               ylabel="",
               xticks=false,
               yticks=false,
               framestyle=:none,
               size=(800, 600),
               dpi=300,
               margin=5Plots.mm,
               right_margin=15Plots.mm,
               clims=(-maximum(abs.(color_matrix)), maximum(abs.(color_matrix))),
               yflip=false)
    

    polynomial = Polynomial(poly_str)

    slp = polynomial.slp

    xs = range(box.x.lower, box.x.upper; length=(round(Int, box.x.upper - box.x.lower) + 1) * 50)
    ys = range(box.y.lower, box.y.upper; length=(round(Int, box.y.upper - box.y.lower) + 1) * 50)
    Z = [ eval_slp(slp, x, y, "") for  y in ys, x in xs ]

    contour!(p, xs, ys, Z, levels=[0], 
             color=:black, linewidth=3, linestyle=:solid,
             label="f(x,y) = 0")

    return p
end
