#=
    evaluate(poly::Polynomial, B::myBox)

    Evaluates the polynomial at the center point of the box.

    # Arguments
    - `poly::Polynomial`: The polynomial
    - `B::myBox`: The box
=#
function evaluate_Taylor(poly::Polynomial, B::myBox, sharing::Bool)
    global derivatives_taylor 
    global total_eval
    global total_points
    global max_x
    global max_y
    global total_degree

    center_pt = (B.x_2, B.y_2)

    if sharing
        if haskey(derivatives_taylor , center_pt)
            return 
        end
    end

    if max_x < 2 || max_y < 2
        point_derivatives = [zeros(Float64, 3) for _ in 1:3]
    else 
        point_derivatives = [zeros(Float64, max_y + 1) for _ in 1:max_x + 1]
    end

    for x_index in range(1, max_x + 1)
        for y_index in range(1, max_y + 1)
            if (x_index-1) + (y_index-1) > total_degree
                continue
            end
            order = monomial_string(x_index-1, y_index-1)
            deriv = eval_slp(poly.slp, center_pt[1], center_pt[2], order) 
            point_derivatives[x_index][y_index] = deriv
            total_eval += 1
        end
    end

    total_points += 1


    derivatives_taylor[center_pt] = point_derivatives
end


@inline quad_eval(c0, c1, c2, t) = (c2*t + c1)*t + c0

@inline function _taylorD(D, ox::Int, oy::Int)::Float64
    xi = ox + 1
    yi = oy + 1
    if xi <= length(D)
        row = D[xi]
        if yi <= length(row)
            return row[yi]
        end
    end
    return 0.0
end

function taylor_interpolation3(poly::Polynomial, B::myBox, degree::Int; sharing::Bool=true)::myInterval
    evaluate_Taylor(poly, B, sharing)

    x_2, y_2 = B.x_2, B.y_2
    rx, ry = B.rx, B.ry

    # This implementation follows the Maple Taylor3 code (square box, centered coords)
    r = rx

    center_pt = (x_2, y_2)
    if !haskey(derivatives_taylor, center_pt)
        return myInterval(0.0, 0.0)
    end
    D = derivatives_taylor[center_pt]

    # coefficients of quadratic Taylor polynomial p in centered form
    c00 = _taylorD(D, 0, 0)
    c10 = _taylorD(D, 1, 0)
    c01 = _taylorD(D, 0, 1)
    c20 = _taylorD(D, 2, 0) / 2
    c11 = _taylorD(D, 1, 1)
    c02 = _taylorD(D, 0, 2) / 2

    # values of p at the corners
    p_aa = c00 - r*(c10+c01) + r^2*(c20+c11+c02)
    p_ab = c00 - r*(c10-c01) + r^2*(c20-c11+c02)
    p_ba = c00 + r*(c10-c01) + r^2*(c20-c11+c02)
    p_bb = c00 + r*(c10+c01) + r^2*(c20+c11+c02)

    alpha = min(p_aa, p_ab, p_ba, p_bb)
    beta  = max(p_aa, p_ab, p_ba, p_bb)

    # analyse quadratic boundary polynomials p1..p4
    eps0 = 1e-14
    for i in 1:4
        if i == 1
            c0 = c00 - r*c01 + r^2*c02
            c1 = c10 - r*c11
            c2 = c20
        elseif i == 2
            c0 = c00 + r*c01 + r^2*c02
            c1 = c10 + r*c11
            c2 = c20
        elseif i == 3
            c0 = c00 - r*c10 + r^2*c20
            c1 = c01 - r*c11
            c2 = c02
        else
            c0 = c00 + r*c10 + r^2*c20
            c1 = c01 + r*c11
            c2 = c02
        end

        if abs(c2) > eps0 && (abs(c1) < 2*abs(c2)*r)
            e = c0 - c1^2/(4*c2)
            if c2 > 0
                alpha = min(alpha, e)
            else
                beta = max(beta, e)
            end
        end
    end

    # analyse p over the interior of B
    Dp = 4*c20*c02 - c11^2
    if Dp > 0
        if (abs(2*c10*c02 - c01*c11) < Dp*r) && (abs(2*c01*c20 - c10*c11) < Dp*r)
            e = c00 - (c10^2*c02 - c10*c01*c11 + c01^2*c20)/Dp
            if c20 > 0
                alpha = min(alpha, e)
            else
                beta = max(beta, e)
            end
        end
    end

    # add remainder bound (k = 3..degree)
    S = 0.0
    for k in 3:degree
        sk = 0.0
        for j in 0:k
            sk += binomial(k, j) * abs(_taylorD(D, k-j, j))
        end
        sk /= float(factorial(big(k)))
        S += sk * r^(k-3)
    end

    return myInterval(alpha - r^3*S, beta + r^3*S)
end
