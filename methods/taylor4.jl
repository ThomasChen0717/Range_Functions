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


@inline cubic_eval(c0,c1,c2,c3,t) = ((c3*t + c2)*t + c1)*t + c0

@inline function update_if_in_interval!(a, b, t, c0,c1,c2,c3, r; eps=1e-12)
    if abs(t) ≤ r + eps
        v = cubic_eval(c0,c1,c2,c3,t)
        a = min(a, v)
        b = max(b, v)   
    end
    return a, b
end

# polynomial coeffs in y, stored low->high: a0 + a1*y + a2*y^2 + ...
@inline function poly_add(a::Vector{Float64}, b::Vector{Float64})
    n = max(length(a), length(b))
    c = zeros(Float64, n)
    for i in 1:n
        c[i] = (i <= length(a) ? a[i] : 0.0) + (i <= length(b) ? b[i] : 0.0)
    end
    return c
end

@inline function poly_mul(a::Vector{Float64}, b::Vector{Float64})
    c = zeros(Float64, length(a) + length(b) - 1)
    for i in 1:length(a), j in 1:length(b)
        c[i + j - 1] += a[i] * b[j]
    end
    return c
end

@inline poly_scale(a::Vector{Float64}, s::Float64) = [s * ai for ai in a]

#=
resultant_quad_quad(A,B,C,D,E,F)

Resultant of two quadratics in x:
f = A x^2 + B x + C
g = D x^2 + E x + F
where A..F are polynomials in y (low->high coeff vectors).
Returns Res(y) as polynomial coeff vector (low->high).
=#
function resultant_quad_quad(A,B,C,D,E,F)
    # Sylvester matrix 4x4 with polynomial entries
    z = [0.0]  # zero polynomial

    # IMPORTANT: matrix of polynomial-vectors, not concatenation
    M = Matrix{Vector{Float64}}(undef, 4, 4)
    M[1,1] = A;  M[1,2] = B;  M[1,3] = C;  M[1,4] = z
    M[2,1] = z;  M[2,2] = A;  M[2,3] = B;  M[2,4] = C
    M[3,1] = D;  M[3,2] = E;  M[3,3] = F;  M[3,4] = z
    M[4,1] = z;  M[4,2] = D;  M[4,3] = E;  M[4,4] = F

    # det via Leibniz (24 perms) since matrix is tiny and entries are polys
    perms = (
        (1,2,3,4,+1),(1,2,4,3,-1),(1,3,2,4,-1),(1,3,4,2,+1),(1,4,2,3,+1),(1,4,3,2,-1),
        (2,1,3,4,-1),(2,1,4,3,+1),(2,3,1,4,+1),(2,3,4,1,-1),(2,4,1,3,-1),(2,4,3,1,+1),
        (3,1,2,4,+1),(3,1,4,2,-1),(3,2,1,4,-1),(3,2,4,1,+1),(3,4,1,2,+1),(3,4,2,1,-1),
        (4,1,2,3,-1),(4,1,3,2,+1),(4,2,1,3,+1),(4,2,3,1,-1),(4,3,1,2,-1),(4,3,2,1,+1),
    )

    res = [0.0]
    for (p1,p2,p3,p4,sgn) in perms
        term = M[1,p1]
        term = poly_mul(term, M[2,p2])
        term = poly_mul(term, M[3,p3])
        term = poly_mul(term, M[4,p4])
        res = poly_add(res, poly_scale(term, float(sgn)))
    end

    # trim tiny trailing coeffs
    while length(res) > 1 && abs(res[end]) < 1e-12
        pop!(res)
    end
    return res
end


function taylor4_boundary!(a, b, c00,c10,c01,c20,c11,c02,c30,c21,c12,c03, r)
    eps0 = 1e-14
    epsD = 1e-14
    epsI = 1e-12

    for i in 1:4
        if i == 1
            c0 = c00 - r*c01 + r^2*c02 - r^3*c03
            c1 = c10 - r*c11 + r^2*c12
            c2 = c20 - r*c21
            c3 = c30
        elseif i == 2
            c0 = c00 + r*c01 + r^2*c02 + r^3*c03
            c1 = c10 + r*c11 + r^2*c12
            c2 = c20 + r*c21
            c3 = c30
        elseif i == 3
            c0 = c00 - r*c10 + r^2*c20 - r^3*c30
            c1 = c01 - r*c11 + r^2*c21
            c2 = c02 - r*c12
            c3 = c03
        else
            c0 = c00 + r*c10 + r^2*c20 + r^3*c30
            c1 = c01 + r*c11 + r^2*c21
            c2 = c02 + r*c12
            c3 = c03
        end

        # endpoints always matter
        a = min(a, cubic_eval(c0,c1,c2,c3,-r), cubic_eval(c0,c1,c2,c3, r))
        b = max(b, cubic_eval(c0,c1,c2,c3,-r), cubic_eval(c0,c1,c2,c3, r))

        if abs(c3) ≤ eps0
            # quadratic: c0 + c1 t + c2 t^2
            if abs(c2) > eps0 && (abs(c1) < 2*abs(c2)*r)
                e = c0 - c1^2/(4*c2)
                if c2 > 0
                    a = min(a, e)
                else
                    b = max(b, e)   
                end
            end
        else
            delta = c2^2 - 3*c1*c3
            if delta > epsD
                L = sign(c3) * (c1 + 3*c3*r^2)
                R = 2*abs(c2)*r
                s_delta = sqrt(delta)

                if L > R
                    if abs(c2) < 3*abs(c3)*r
                        xm = -(c2 - s_delta)/(3*c3)
                        xp = -(c2 + s_delta)/(3*c3)
                        a, b = update_if_in_interval!(a, b, xm, c0,c1,c2,c3, r; eps=epsI)
                        a, b = update_if_in_interval!(a, b, xp, c0,c1,c2,c3, r; eps=epsI)
                    end
                elseif L > -R
                    if c2 > 0
                        xm = -(c2 - s_delta)/(3*c3)
                        a, b = update_if_in_interval!(a, b, xm, c0,c1,c2,c3, r; eps=epsI)
                    elseif c2 < 0
                        xp = -(c2 + s_delta)/(3*c3)
                        a, b = update_if_in_interval!(a, b, xp, c0,c1,c2,c3, r; eps=epsI)
                    end
                end
            end
        end
    end

    return a, b
end


function taylor_interpolation4(poly::Polynomial, B::myBox, degree::Int; sharing::Bool=true)::myInterval
    evaluate_Taylor(poly, B, sharing)

    x_2, y_2 = B.x_2, B.y_2
    rx, ry = B.rx, B.ry

    r = rx

    center_pt = (x_2, y_2)
    if !haskey(derivatives_taylor, center_pt)
        return myInterval(0.0, 0.0)
    end
    D = derivatives_taylor[center_pt]

    # cubic Taylor coefficients in centered coords (x,y are offsets from center)
    c00 = D[1][1]
    c10 = D[2][1]
    c01 = D[1][2]
    c20 = D[3][1] / 2
    c11 = D[2][2]
    c02 = D[1][3] / 2
    c30 = D[4][1] / 6
    c21 = D[3][2] / 2
    c12 = D[2][3] / 2
    c03 = D[1][4] / 6

    # corners of [-r,r]^2 in centered coordinates
    p_aa = c00 - r*(c10+c01) + r^2*(c20+c11+c02) - r^3*(c30+c21+c12+c03)
    p_ab = c00 - r*(c10-c01) + r^2*(c20-c11+c02) - r^3*(c30-c21+c12-c03)
    p_ba = c00 + r*(c10-c01) + r^2*(c20-c11+c02) + r^3*(c30-c21+c12-c03)
    p_bb = c00 + r*(c10+c01) + r^2*(c20+c11+c02) + r^3*(c30+c21+c12+c03)

    a = min(p_aa,p_ab,p_ba,p_bb)
    b = max(p_aa,p_ab,p_ba,p_bb)

    # boundary cubics (p1..p4)
    a, b = taylor4_boundary!(a, b, c00,c10,c01,c20,c11,c02,c30,c21,c12,c03, r)

    # interior stationary points
    a, b = interior!(a, b, c00,c10,c01,c20,c11,c02,c30,c21,c12,c03, r)

    # remainder bound (k=4..degree)
    S = 0.0
    for k in 4:degree
        sk = 0.0
        for j in 0:k
            dk = D[(k-j)+1][j+1]
            sk += binomial(k,j) * abs(dk)
        end
        sk /= float(factorial(big(k)))
        S += sk * r^(k-4)
    end

    return myInterval(a - r^4*S, b + r^4*S)
end


