#=
Lagrange3.jl

Translation of the provided Maple `Lagrange3` routine into Julia.

Primary functions:
- `trinomial(n,k)` : integer trinomial coefficient (matches the Maple recursion)
- `RangeBi2(c00,c10,c01,c20,c11,c02,c21,c12,c22,r)` : computes range of a biquadratic polynomial in centered form over box radius `r` (returns (alpha,beta))
- `mixed_partial(f, x0, y0, ox, oy)` : computes mixed partial derivative ∂^(ox+oy) f / ∂x^ox ∂y^oy at (x0,y0) using ForwardDiff
- `Lagrange3(f, B; max_k=1)` : main routine. `f` is `f(x,y)`, `B` is tuple `(mx,my,r)`. `max_k` is the number of remainder terms to include (Maple used floor(degree/3)); default 1.
=# 
const TRIN_TABLE = Ref(zeros(Int, 1, 1))
const TRIN_MAX   = Ref(0)
const LAGRANGE_ORDER_CACHE = Dict{Tuple{Int, Int, Int}, Vector{Tuple{Int, Int, String}}}()

struct UniformLagrangeGridCache
    xmin::Float64
    ymin::Float64
    step::Float64
    r::Float64
    order_index::Matrix{Int}
    values::Vector{Matrix{Float64}}
    remainder_terms::Vector{Tuple{Int, Float64}}
end

function trinomial_table!(nmax::Int)
    nmax <= TRIN_MAX[] && return

    T = zeros(Int, nmax+1, nmax+1)
    for n in 0:nmax
        T[n+1, 0+1] = 1
        T[n+1, n+1] = 1
        for k in 1:(n-1)
            t1 = T[n,   k]     # T(n-1,k-1)
            t2 = T[n,   k+1]   # T(n-1,k)
            t3 = n >= 2 ? T[n-1, k] : 0  # T(n-2,k-1)
            T[n+1, k+1] = t1 + t2 + t3
        end
    end

    TRIN_TABLE[] = T
    TRIN_MAX[]   = nmax
end

function trinomial(n::Int, k::Int)
    (k < 0 || k > n) && return 0
    n > TRIN_MAX[] && trinomial_table!(n)
    @inbounds return TRIN_TABLE[][n+1, k+1]
end

function lagrange_orders(max_x::Int, max_y::Int, total_degree::Int)::Vector{Tuple{Int, Int, String}}
    key = (max_x, max_y, total_degree)
    return get!(LAGRANGE_ORDER_CACHE, key) do
        max_x_multiples = div(max_x, 3) + 1
        max_y_multiples = div(max_y, 3) + 1
        orders = Tuple{Int, Int, String}[]

        for x_index in 1:max_x_multiples
            for y_index in 1:max_y_multiples
                x_order = (x_index - 1) * 3
                y_order = (y_index - 1) * 3
                if x_order + y_order > total_degree
                    continue
                end
                push!(orders, (x_index, y_index, monomial_string(x_order, y_order)))
            end
        end

        orders
    end
end

function build_uniform_lagrange_grid_cache(poly::Polynomial, boxes::Vector{myBox}, degree::Int)
    isempty(boxes) && return nothing

    first_box = boxes[1]
    tol = 1e-12
    if !isapprox(first_box.rx, first_box.ry; atol=tol, rtol=tol)
        return nothing
    end

    for box in boxes
        if !isapprox(box.rx, first_box.rx; atol=tol, rtol=tol) || !isapprox(box.ry, first_box.ry; atol=tol, rtol=tol)
            return nothing
        end
    end

    xmin = minimum(box.x_1 for box in boxes)
    xmax = maximum(box.x_3 for box in boxes)
    ymin = minimum(box.y_1 for box in boxes)
    ymax = maximum(box.y_3 for box in boxes)
    step = first_box.rx

    nx = round(Int, (xmax - xmin) / step) + 1
    ny = round(Int, (ymax - ymin) / step) + 1

    if nx <= 0 || ny <= 0
        return nothing
    end

    if !isapprox((nx - 1) * step, xmax - xmin; atol=1e-9, rtol=1e-9) || !isapprox((ny - 1) * step, ymax - ymin; atol=1e-9, rtol=1e-9)
        return nothing
    end

    max_x_multiples = div(max_x, 3) + 1
    max_y_multiples = div(max_y, 3) + 1
    order_index = zeros(Int, max_x_multiples, max_y_multiples)
    orders = lagrange_orders(max_x, max_y, total_degree)

    values = Matrix{Float64}[]
    sizehint!(values, length(orders))
    for (x_index, y_index, _) in orders
        push!(values, Matrix{Float64}(undef, nx, ny))
        order_index[x_index, y_index] = length(values)
    end

    vars = Dict{Symbol,Union{Float64,myInterval}}(:x => 0.0, :y => 0.0)
    slp = poly.slp
    for ix in 1:nx
        vars[:x] = xmin + (ix - 1) * step
        for iy in 1:ny
            vars[:y] = ymin + (iy - 1) * step
            for (x_index, y_index, order) in orders
                idx = order_index[x_index, y_index]
                @inbounds values[idx][ix, iy] = evaluate_slp_range(slp, order, vars)
            end
        end
    end

    max_k = degree ÷ 3
    trinomial_table!(max_k)
    omega = sqrt(3) / 27 * first_box.rx^3
    remainder_terms = Tuple{Int, Float64}[]
    for k in 1:max_k
        for j in 0:k
            idx = order_index[k - j + 1, j + 1]
            if idx == 0
                continue
            end
            push!(remainder_terms, (idx, omega^k * trinomial(k, j)))
        end
    end

    return UniformLagrangeGridCache(xmin, ymin, step, first_box.rx, order_index, values, remainder_terms)
end

@inline function lattice_box_origin(cache::UniformLagrangeGridCache, B::myBox)
    ix = round(Int, (B.x_1 - cache.xmin) / cache.step) + 1
    iy = round(Int, (B.y_1 - cache.ymin) / cache.step) + 1
    return ix, iy
end

function Lagrange3(cache::UniformLagrangeGridCache, B::myBox)::myInterval
    ix, iy = lattice_box_origin(cache, B)
    r = cache.r

    base = cache.values[cache.order_index[1, 1]]
    @inbounds begin
        f11 = base[ix, iy]
        f12 = base[ix, iy + 1]
        f13 = base[ix, iy + 2]
        f21 = base[ix + 1, iy]
        f22 = base[ix + 1, iy + 1]
        f23 = base[ix + 1, iy + 2]
        f31 = base[ix + 2, iy]
        f32 = base[ix + 2, iy + 1]
        f33 = base[ix + 2, iy + 2]

        c00 = f22
        c10 = (f32 - f12) / (2*r)
        c01 = (f23 - f21) / (2*r)
        c20 = (f32 - 2*f22 + f12) / (2*r^2)
        c11 = (f33 - f13 - f31 + f11) / (4*r^2)
        c02 = (f23 - 2*f22 + f21) / (2*r^2)
        c21 = (f33 - 2*f23 + f13 - f31 + 2*f21 - f11) / (4*r^3)
        c12 = (f33 - 2*f32 + f31 - f13 + 2*f12 - f11) / (4*r^3)
        c22 = (f33 - 2*f23 + f13 - 2*f32 + 4*f22 - 2*f12 + f31 - 2*f21 + f11) / (4*r^4)

        U0 = RangeBi2(c00, c10, c01, c20, c11, c02, c21, c12, c22, r)

        U3_alpha = 0.0
        U3_beta = 0.0
        for (idx, coeff) in cache.remainder_terms
            mat = cache.values[idx]

            g11 = mat[ix, iy]
            g12 = mat[ix, iy + 1]
            g13 = mat[ix, iy + 2]
            g21 = mat[ix + 1, iy]
            g22 = mat[ix + 1, iy + 1]
            g23 = mat[ix + 1, iy + 2]
            g31 = mat[ix + 2, iy]
            g32 = mat[ix + 2, iy + 1]
            g33 = mat[ix + 2, iy + 2]

            cc00 = g22
            cc10 = (g32 - g12) / (2*r)
            cc01 = (g23 - g21) / (2*r)
            cc20 = (g32 - 2*g22 + g12) / (2*r^2)
            cc11 = (g33 - g13 - g31 + g11) / (4*r^2)
            cc02 = (g23 - 2*g22 + g21) / (2*r^2)
            cc21 = (g33 - 2*g23 + g13 - g31 + 2*g21 - g11) / (4*r^3)
            cc12 = (g33 - 2*g32 + g31 - g13 + 2*g12 - g11) / (4*r^3)
            cc22 = (g33 - 2*g23 + g13 - 2*g32 + 4*g22 - 2*g12 + g31 - 2*g21 + g11) / (4*r^4)

            u_alpha, u_beta = RangeBi2(cc00, cc10, cc01, cc20, cc11, cc02, cc21, cc12, cc22, r)
            term = coeff * max(abs(u_alpha), abs(u_beta))
            U3_alpha += term
            U3_beta += term
        end

        return myInterval(U0[1] - U3_alpha, U0[2] + U3_beta)
    end
end

#= 
    evaluate(poly::Polynomial, B::myBox)

    Evaluates the polynomial and computes derivatives

    # Arguments
    - `poly`: The polynomial to evaluate
    - `B`: The box region where the evaluation is performed
=#
function evaluate(poly::Polynomial, B::myBox, sharing::Bool)
    slp = poly.slp


    global derivatives
    global total_points
    global total_eval
    global max_x
    global max_y
    global total_degree

    max_x_multiples = div(max_x, 3) + 1
    max_y_multiples = div(max_y, 3) + 1
    orders = lagrange_orders(max_x, max_y, total_degree)
    vars = Dict{Symbol,Union{Float64,myInterval}}(:x => 0.0, :y => 0.0)

    for i in range(1, 3)
        for j in range(1,3)
            pt = get_point(B, i, j)

            if sharing
                if haskey(derivatives, pt)
                    continue 
                end
            end

            point_derivatives = [zeros(Float64, max_y_multiples) for _ in 1:max_x_multiples]
            vars[:x] = pt[1]
            vars[:y] = pt[2]

            for (x_index, y_index, order) in orders
                    deriv = evaluate_slp_range(slp, order, vars)
                    point_derivatives[x_index][y_index] = deriv
                    total_eval += 1
            end

            total_points += 1
            derivatives[pt] = point_derivatives
        end
    end
end

@inline function get_cached_deriv(pt_derivs, ox::Int, oy::Int)
    # only multiples of 3 are cached
    (ox < 0 || oy < 0 || ox % 3 != 0 || oy % 3 != 0) && return 0.0

    xi = ox ÷ 3 + 1
    yi = oy ÷ 3 + 1

    (xi > length(pt_derivs)) && return 0.0
    (yi > length(pt_derivs[xi])) && return 0.0

    @inbounds return pt_derivs[xi][yi]
end

@inline function cached_mixed_partial(pt::Tuple{Float64,Float64}, ox::Int, oy::Int)
    global derivatives
    haskey(derivatives, pt) || error("No cached derivatives for point $pt. Did you call evaluate(poly, B, sharing=true) first?")
    return get_cached_deriv(derivatives[pt], ox, oy)
end


function RangeBi2(c00,c10,c01,c20,c11,c02,c21,c12,c22,r)
    # Translate Maple RangeBi2 logic to Julia
    q_aa = c00 - r*(c10 + c01) + r^2*(c20 + c11 + c02)
    q_ab = c00 - r*(c10 - c01) + r^2*(c20 - c11 + c02)
    q_ba = c00 + r*(c10 - c01) + r^2*(c20 - c11 + c02)
    q_bb = c00 + r*(c10 + c01) + r^2*(c20 + c11 + c02)

    alpha_q = minimum((q_aa,q_ab,q_ba,q_bb))
    beta_q  = maximum((q_aa,q_ab,q_ba,q_bb))

    # analyse quadratic boundary polynomials q1..q4
    for i in 1:4
        if i == 1
            c0 = c00 - r*c01 + r^2*c02; c1 = c10 - r*c11; c2 = c20
        elseif i == 2
            c0 = c00 + r*c01 + r^2*c02; c1 = c10 + r*c11; c2 = c20
        elseif i == 3
            c0 = c00 - r*c10 + r^2*c20; c1 = c01 - r*c11; c2 = c02
        else
            c0 = c00 + r*c10 + r^2*c20; c1 = c01 + r*c11; c2 = c02
        end

        if abs(c1) < 2*abs(c2)*r && c2 != 0
            e = c0 - c1^2 / (4*c2)
            if c2 > 0
                alpha_q = min(alpha_q, e)
            else
                beta_q  = max(beta_q, e)
            end
        end
    end

    # analyse q over interior
    Dq = 4*c20*c02 - c11^2
    if Dq > 0
        cond1 = abs(2*c10*c02 - c01*c11) < Dq*r
        cond2 = abs(2*c01*c20 - c10*c11) < Dq*r
        if cond1 && cond2
            e = c00 - (c10^2*c02 - c10*c01*c11 + c01^2*c20)/Dq
            if c20 > 0
                alpha_q = min(alpha_q, e)
            else
                beta_q  = max(beta_q, e)
            end
        end
    end

    # remainder polynomial r values at corners
    r_aa = r^4*c22 - r^3*(c21 + c12)
    r_ab = r^4*c22 - r^3*(c21 - c12)
    r_ba = r^4*c22 + r^3*(c21 - c12)
    r_bb = r^4*c22 + r^3*(c21 + c12)

    alpha_r = minimum((r_aa,r_ab,r_ba,r_bb))
    beta_r  = maximum((r_aa,r_ab,r_ba,r_bb))

    # analyse quadratic boundaries for r
    for i in 1:4
        if i == 1
            c0 = 0; c1 = c12*r^2; c2 = c22*r^2 - c21*r
        elseif i == 2
            c0 = 0; c1 = c12*r^2; c2 = c22*r^2 + c21*r
        elseif i == 3
            c0 = 0; c1 = c21*r^2; c2 = c22*r^2 - c12*r
        else
            c0 = 0; c1 = c21*r^2; c2 = c22*r^2 + c12*r
        end

        if abs(c1) < 2*abs(c2)*r && c2 != 0
            e = c0 - c1^2 / (4*c2)
            if c2 > 0
                alpha_r = min(alpha_r, e)
            else
                beta_r  = max(beta_r, e)
            end
        end
    end

    return (alpha_q + alpha_r, beta_q + beta_r)
end

function Lagrange3(poly::Polynomial, B::myBox, degree::Int; sharing::Bool=true)::myInterval
    # B = (mx, my, r)
    global max_x, max_y

    mx, my = B.x_2, B.y_2
    rx, ry = B.rx, B.ry

    r = rx

    slp = poly.slp

    evaluate(poly, B, sharing)

    max_k = degree ÷ 3
    trinomial_table!(max_k)

    pt11 = get_point(B, 1, 1)
    pt12 = get_point(B, 1, 2)
    pt13 = get_point(B, 1, 3)
    pt21 = get_point(B, 2, 1)
    pt22 = get_point(B, 2, 2)
    pt23 = get_point(B, 2, 3)
    pt31 = get_point(B, 3, 1)
    pt32 = get_point(B, 3, 2)
    pt33 = get_point(B, 3, 3)

    @assert haskey(derivatives, pt11) "Missing cached derivatives at $pt11. Did you call evaluate(poly, B, sharing=true) first?"
    @assert haskey(derivatives, pt12) "Missing cached derivatives at $pt12. Did you call evaluate(poly, B, sharing=true) first?"
    @assert haskey(derivatives, pt13) "Missing cached derivatives at $pt13. Did you call evaluate(poly, B, sharing=true) first?"
    @assert haskey(derivatives, pt21) "Missing cached derivatives at $pt21. Did you call evaluate(poly, B, sharing=true) first?"
    @assert haskey(derivatives, pt22) "Missing cached derivatives at $pt22. Did you call evaluate(poly, B, sharing=true) first?"
    @assert haskey(derivatives, pt23) "Missing cached derivatives at $pt23. Did you call evaluate(poly, B, sharing=true) first?"
    @assert haskey(derivatives, pt31) "Missing cached derivatives at $pt31. Did you call evaluate(poly, B, sharing=true) first?"
    @assert haskey(derivatives, pt32) "Missing cached derivatives at $pt32. Did you call evaluate(poly, B, sharing=true) first?"
    @assert haskey(derivatives, pt33) "Missing cached derivatives at $pt33. Did you call evaluate(poly, B, sharing=true) first?"

    d11 = derivatives[pt11]
    d12 = derivatives[pt12]
    d13 = derivatives[pt13]
    d21 = derivatives[pt21]
    d22 = derivatives[pt22]
    d23 = derivatives[pt23]
    d31 = derivatives[pt31]
    d32 = derivatives[pt32]
    d33 = derivatives[pt33]

    f11 = d11[1][1]
    f12 = d12[1][1]
    f13 = d13[1][1]
    f21 = d21[1][1]
    f22 = d22[1][1]
    f23 = d23[1][1]
    f31 = d31[1][1]
    f32 = d32[1][1]
    f33 = d33[1][1]

    c00 = f22
    c10 = (f32 - f12) / (2*r)
    c01 = (f23 - f21) / (2*r)
    c20 = (f32 - 2*f22 + f12) / (2*r^2)
    c11 = (f33 - f13 - f31 + f11) / (4*r^2)
    c02 = (f23 - 2*f22 + f21) / (2*r^2)
    c21 = (f33 - 2*f23 + f13 - f31 + 2*f21 - f11) / (4*r^3)
    c12 = (f33 - 2*f32 + f31 - f13 + 2*f12 - f11) / (4*r^3)
    c22 = (f33 - 2*f23 + f13 - 2*f32 + 4*f22 - 2*f12 + f31 - 2*f21 + f11) / (4*r^4)

    U0 = RangeBi2(c00,c10,c01,c20,c11,c02,c21,c12,c22,r)

    Omega = sqrt(3)/27 * r^3
    U3_alpha = 0.0
    U3_beta = 0.0

    # sum remainder terms up to max_k (Maple used n=floor(degree/3)).
    for k in 1:max_k
        for j in 0:k
            ox = 3*(k-j)
            oy = 3*j

            g11 = get_cached_deriv(d11, ox, oy)
            g12 = get_cached_deriv(d12, ox, oy)
            g13 = get_cached_deriv(d13, ox, oy)
            g21 = get_cached_deriv(d21, ox, oy)
            g22 = get_cached_deriv(d22, ox, oy)
            g23 = get_cached_deriv(d23, ox, oy)
            g31 = get_cached_deriv(d31, ox, oy)
            g32 = get_cached_deriv(d32, ox, oy)
            g33 = get_cached_deriv(d33, ox, oy)

            cc00 = g22
            cc10 = (g32 - g12) / (2*r)
            cc01 = (g23 - g21) / (2*r)
            cc20 = (g32 - 2*g22 + g12) / (2*r^2)
            cc11 = (g33 - g13 - g31 + g11) / (4*r^2)
            cc02 = (g23 - 2*g22 + g21) / (2*r^2)
            cc21 = (g33 - 2*g23 + g13 - g31 + 2*g21 - g11) / (4*r^3)
            cc12 = (g33 - 2*g32 + g31 - g13 + 2*g12 - g11) / (4*r^3)
            cc22 = (g33 - 2*g23 + g13 - 2*g32 + 4*g22 - 2*g12 + g31 - 2*g21 + g11) / (4*r^4)

            u_alpha, u_beta = RangeBi2(cc00,cc10,cc01,cc20,cc11,cc02,cc21,cc12,cc22,r)
            term = Omega^k * trinomial(k, j) * max(abs(u_alpha), abs(u_beta))
            U3_alpha += term
            U3_beta  += term
        end
    end

    return myInterval(U0[1] - U3_alpha, U0[2] + U3_beta)
end

