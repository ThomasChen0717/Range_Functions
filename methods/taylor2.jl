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

function _taylorD2(D, ox::Int, oy::Int)::Float64
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

function taylor_interpolation2(poly::Polynomial, B::myBox, degree::Int; sharing::Bool=true)::myInterval
    # Implements Maple Taylor2
    # Returns [f(mx,my) - r*S, f(mx,my) + r*S]
    # where S = sum_{k=1..d} s[k]*r^(k-1),
    # s[k] = (1/k!) * sum_{j=0..k} binomial(k,j)*abs(df[k-j,j](mx,my)).

    evaluate_Taylor(poly, B, sharing)

    x_2, y_2 = B.x_2, B.y_2
    rx, ry = B.rx, B.ry

    # Follow existing Taylor code: assume square, use rx
    r = rx

    center_pt = (x_2, y_2)
    if !haskey(derivatives_taylor, center_pt)
        return myInterval(0.0, 0.0)
    end
    D = derivatives_taylor[center_pt]

    s0 = _taylorD2(D, 0, 0)

    S = 0.0
    for k in 1:degree
        sk = 0.0
        for j in 0:k
            sk += binomial(k, j) * abs(_taylorD2(D, k-j, j))
        end
        sk /= float(factorial(big(k)))
        S += sk * r^(k-1)
    end

    return myInterval(s0 - r*S, s0 + r*S)
end
