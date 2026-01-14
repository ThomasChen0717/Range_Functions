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

    for i in range(1, 3)
        for j in range(1,3)
            pt = get_point(B, i, j)

            if sharing
                if haskey(derivatives, pt)
                    continue 
                end
            end

            max_x_multiples = div(max_x, 3) + 1
            max_y_multiples = div(max_y, 3) + 1  

            point_derivatives = [zeros(Float64, max_y_multiples+ 1) for _ in 1:max_x_multiples + 1]

            for x_index in range(1, max_x_multiples)
                for y_index in range(1, max_y_multiples)
                    x_order = (x_index - 1) * 3
                    y_order = (y_index - 1) * 3
                    if x_order + y_order > total_degree
                        continue
                    end
                    order = monomial_string(x_order, y_order)
                    deriv = eval_slp(slp, pt[1], pt[2], order) 
                    point_derivatives[x_index][y_index] = deriv
                    total_eval += 1
                end
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

    # build F grid values F[0..2,0..2] where index i corresponds to i-1 in Maple
    F = zeros(Float64, 3, 3)

    for i in 1:3, j in 1:3
        pt = get_point(B, i, j)  

        @assert haskey(derivatives, pt) "Missing cached derivatives at $pt. Did you call evaluate(poly, B, sharing=true) first?"
        
        F[i, j] = derivatives[pt][1][1]
    end

    c00 = F[2,2]
    c10 = (F[3,2] - F[1,2]) / (2*r)
    c01 = (F[2,3] - F[2,1]) / (2*r)
    c20 = (F[3,2] - 2*F[2,2] + F[1,2]) / (2*r^2)
    c11 = (F[3,3] - F[1,3] - F[3,1] + F[1,1]) / (4*r^2)
    c02 = (F[2,3] - 2*F[2,2] + F[2,1]) / (2*r^2)
    c21 = (F[3,3] - 2*F[2,3] + F[1,3] - F[3,1] + 2*F[2,1] - F[1,1]) / (4*r^3)
    c12 = (F[3,3] - 2*F[3,2] + F[3,1] - F[1,3] + 2*F[1,2] - F[1,1]) / (4*r^3)
    c22 = (F[3,3] - 2*F[2,3] + F[1,3] - 2*F[3,2] + 4*F[2,2] - 2*F[1,2] + F[3,1] - 2*F[2,1] + F[1,1]) / (4*r^4)

    U0 = RangeBi2(c00,c10,c01,c20,c11,c02,c21,c12,c22,r)

    Omega = sqrt(3)/27 * r^3
    U3_alpha = 0.0
    U3_beta = 0.0

    # sum remainder terms up to max_k (Maple used n=floor(degree/3)).
    for k in 1:max_k
        for j in 0:k
            ox = 3*(k-j)
            oy = 3*j

            # compute 3(k-j), 3j derivatives at the 3x3 grid
            G = zeros(Float64, 3, 3)

            for ii in 1:3, jj in 1:3
                pt = get_point(B, ii, jj)  # Tuple{Float64,Float64}

                @assert haskey(derivatives, pt) "Missing cached derivatives at $pt. Did you call evaluate(poly, B, sharing=true) first?"

                # ox, oy are multiples of 3; cached at [ox/3+1][oy/3+1]
                G[ii, jj] = get_cached_deriv(derivatives[pt], ox, oy)
            end

            # coefficients of biquadratic Lagrange polynomial for this derivative array
            cc00 = G[2,2]
            cc10 = (G[3,2] - G[1,2]) / (2*r)
            cc01 = (G[2,3] - G[2,1]) / (2*r)
            cc20 = (G[3,2] - 2*G[2,2] + G[1,2]) / (2*r^2)
            cc11 = (G[3,3] - G[1,3] - G[3,1] + G[1,1]) / (4*r^2)
            cc02 = (G[2,3] - 2*G[2,2] + G[2,1]) / (2*r^2)
            cc21 = (G[3,3] - 2*G[2,3] + G[1,3] - G[3,1] + 2*G[2,1] - G[1,1]) / (4*r^3)
            cc12 = (G[3,3] - 2*G[3,2] + G[3,1] - G[1,3] + 2*G[1,2] - G[1,1]) / (4*r^3)
            cc22 = (G[3,3] - 2*G[2,3] + G[1,3] - 2*G[3,2] + 4*G[2,2] - 2*G[1,2] + G[3,1] - 2*G[2,1] + G[1,1]) / (4*r^4)

            u_alpha, u_beta = RangeBi2(cc00,cc10,cc01,cc20,cc11,cc02,cc21,cc12,cc22,r)
            term = Omega^k * trinomial(k, j) * max(abs(u_alpha), abs(u_beta))
            U3_alpha += term
            U3_beta  += term
        end
    end

    return myInterval(U0[1] - U3_alpha, U0[2] + U3_beta)
end

# B = (0.0, 0.0, 1.0)
# poly_str = "x^4 + 2x^2*y^2 + y^4"

# poly = Polynomial(poly_str)
# compute_third_derivatives_2D!(poly)

# total_degree = get_total_degree(poly)
# max_k = total_degree ÷ 3

# print("Range = ", Lagrange3(poly, B, max_k=max_k))
