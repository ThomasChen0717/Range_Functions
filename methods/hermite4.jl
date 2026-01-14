const TRIN_TABLE = Ref(zeros(Int, 1, 1))
const TRIN_MAX   = Ref(0)
const PERMS5 = Ref{Vector{NTuple{6,Int}}}(NTuple{6,Int}[])

@inline normpt(x::Float64, y::Float64) = (round(x, digits=15), round(y, digits=15))

# @inline function get_cached_deriv(pt_derivs::Vector{Vector{Float64}},
#                                   a::Int, b::Int,
#                                   max_x::Int, max_y::Int,
#                                   degree::Int)
#     # outside polynomial degree constraint
#     (a < 0 || b < 0 || a + b > degree) && return 0.0
#     # outside what you cached
#     (a > max_x || b > max_y) && return 0.0
#     return pt_derivs[a+1][b+1]
# end

@inline function get_cached_deriv(pt_derivs, a, b, max_x, max_y, degree)
    if(a < 0 || b < 0 || a + b > degree) 
        return 0.0
    end
    if(a > max_x || b > max_y) 
        return 0.0
    end
    v = pt_derivs[a+1][b+1]
    @assert !isnan(v) "Missing cached derivative (a,b)=($a,$b)"
    return v
end

function evaluate_hermite(poly::Polynomial, B::myBox, sharing::Bool)
    slp = poly.slp

    global derivatives_hermite
    global total_degree
    global max_x
    global max_y

    # We will store all (a,b) with 0 ≤ a,b ≤ total_degree
    Nx = max_x + 1   # a = 0..total_degree
    Ny = max_y + 1   # b = 0..total_degree

    n = fld(total_degree, 4)  

    for i in (1, 3)
        for j in (1, 3)
            pt = get_point(B, i, j)

            if sharing && haskey(derivatives_hermite, pt)
                continue
            end

            # Full derivative table at this corner:
            # point_derivatives[ia][ib] corresponds to derivative (a=ia-1, b=ib-1)
            # point_derivatives = [zeros(Float64, Ny) for _ in 1:Nx]
            point_derivatives = [fill(NaN, Ny) for _ in 1:Nx]

            # remainder derivatives needed by hermite4:
            for k in 0:n
                for jj in 0:k
                    ub = k - jj
                    vb = jj

                    xi0 = 4*ub
                    xi1 = xi0 + 1
                    yi0 = 4*vb
                    yi1 = yi0 + 1
                    

                    for (a,b) in ((xi0,yi0), (xi1,yi0), (xi0,yi1), (xi1,yi1))
                        if a ≤ max_x && b ≤ max_y && (a + b) ≤ total_degree
                            point_derivatives[a+1][b+1] =
                                eval_slp(slp, pt[1], pt[2], monomial_string(a, b))
                        end
                    end
                end
            end

            derivatives_hermite[pt] = point_derivatives
        end
    end
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


function q_val(xx, yy, c00,c10,c01,c20,c11,c02,c30,c21,c12,c03)
    return c00 +
           c10*xx + c01*yy +
           c20*xx^2 + c11*xx*yy + c02*yy^2 +
           c30*xx^3 + c21*xx^2*yy + c12*xx*yy^2 + c03*yy^3
end

function r_val(xx, yy, c31,c22,c13,c32,c23,c33)
    return c31*xx^3*yy   +
           c22*xx^2*yy^2 +
           c13*xx*yy^3   +
           c32*xx^3*yy^2 +
           c23*xx^2*yy^3 +
           c33*xx^3*yy^3
end

function _gen_perms5!()
    # stores (p1,p2,p3,p4,p5,sgn) with sgn ∈ {+1,-1}
    perms = NTuple{6,Int}[]
    p = collect(1:5)

    function sign_of_perm(v)
        invs = 0
        for i in 1:5, j in i+1:5
            invs += (v[i] > v[j])
        end
        return isodd(invs) ? -1 : +1
    end

    function rec!(i)
        if i == 6
            s = sign_of_perm(p)
            push!(perms, (p[1],p[2],p[3],p[4],p[5], s))
            return
        end
        for j in i:5
            p[i], p[j] = p[j], p[i]
            rec!(i+1)
            p[i], p[j] = p[j], p[i]
        end
    end

    rec!(1)
    PERMS5[] = perms
end

function resultant_quad_cubic(A2,A1,A0, B3,B2,B1,B0)
    isempty(PERMS5[]) && _gen_perms5!()
    z = [0.0]

    # 5×5 Sylvester matrix with polynomial entries
    M = Matrix{Vector{Float64}}(undef, 5, 5)
    # rows for g shifted (deg 3 => 2 shifts)
    M[1,1]=B3; M[1,2]=B2; M[1,3]=B1; M[1,4]=B0; M[1,5]=z
    M[2,1]=z;  M[2,2]=B3; M[2,3]=B2; M[2,4]=B1; M[2,5]=B0

    # rows for f shifted (deg 2 => 3 shifts)
    M[3,1]=A2; M[3,2]=A1; M[3,3]=A0; M[3,4]=z;  M[3,5]=z
    M[4,1]=z;  M[4,2]=A2; M[4,3]=A1; M[4,4]=A0; M[4,5]=z
    M[5,1]=z;  M[5,2]=z;  M[5,3]=A2; M[5,4]=A1; M[5,5]=A0

    # determinant by permutations (120 terms), polynomial arithmetic
    res = [0.0]
    for (p1,p2,p3,p4,p5,sgn) in PERMS5[]
        term = M[1,p1]
        term = poly_mul(term, M[2,p2])
        term = poly_mul(term, M[3,p3])
        term = poly_mul(term, M[4,p4])
        term = poly_mul(term, M[5,p5])
        res  = poly_add(res, poly_scale(term, float(sgn)))
    end

    while length(res) > 1 && abs(res[end]) < 1e-12
        pop!(res)
    end
    return res
end



function update_alpha_beta_r_interior!(alpha_r, beta_r,
    c31,c22,c13,c32,c23,c33, r)

    # rx = x^2*(3c31 y + 3c32 y^2 + 3c33 y^3) + x*(2c22 y^2 + 2c23 y^3) + (c13 y^3)
    A2 = [0.0, 3c31, 3c32, 3c33]         # 0 + 3c31 y + 3c32 y^2 + 3c33 y^3
    A1 = [0.0, 0.0, 2c22, 2c23]          # 2c22 y^2 + 2c23 y^3
    A0 = [0.0, 0.0, 0.0, c13]            # c13 y^3

    # ry = x^3*(c31 + 2c32 y + 3c33 y^2) + x^2*(2c22 y + 3c23 y^2) + x*(3c13 y^2) + 0
    B3 = [c31, 2c32, 3c33]
    B2 = [0.0, 2c22, 3c23]
    B1 = [0.0, 0.0, 3c13]
    B0 = [0.0]

    # Degeneracy guard (keeps you safe; you can refine later)
    if maximum(abs.(A2)) < 1e-14 || maximum(abs.(B3)) < 1e-14
        return alpha_r, beta_r
    end

    Res = resultant_quad_cubic(A2,A1,A0, B3,B2,B1,B0)
    ys  = real_roots_poly(Res)

    for y in ys
        abs(y) > r + 1e-12 && continue

        # Solve rx(x,y)=0 as quadratic in x
        qa = (3c31*y + 3c32*y^2 + 3c33*y^3)
        qb = (2c22*y^2 + 2c23*y^3)
        qc = (c13*y^3)

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
            abs(x) > r + 1e-12 && continue

            # Check ry ≈ 0
            ryv = (c31 + 2c32*y + 3c33*y^2)*x^3 +
                  (2c22*y + 3c23*y^2)*x^2 +
                  (3c13*y^2)*x

            abs(ryv) > 1e-8 && continue

            v = r_val(x,y, c31,c22,c13,c32,c23,c33)
            alpha_r = min(alpha_r, v)
            beta_r  = max(beta_r,  v)
        end
    end

    return alpha_r, beta_r
end




#=
    range_bi3(c00,c10,c01,c20,c11,c02,c30,c21,c12,c03,c31,c22,c13,c32,c23,c33,radius)

Range of a bicubic polynomial `p` with coefficients c_ij in centered form
over the box [-r, r] x [-r, r], using the split form p = q + r.
Returns a tuple (α, β).
=#
function range_bi3(
    c00,c10,c01,c20,c11,c02,c30,c21,c12,c03,c31,c22,c13,c32,c23,c33,
    r::Real,
) 
    # split p = q+r, where q is the cubic part and r is the remainder

    # values of q at the corners of B
    q_aa=c00 - r*(c10+c01) + r^2*(c20+c11+c02) - r^3*(c30+c21+c12+c03)
    q_ab=c00 - r*(c10-c01) + r^2*(c20-c11+c02) - r^3*(c30-c21+c12-c03)
    q_ba=c00 + r*(c10-c01) + r^2*(c20-c11+c02) + r^3*(c30-c21+c12-c03)
    q_bb=c00 + r*(c10+c01) + r^2*(c20+c11+c02) + r^3*(c30+c21+c12+c03)

    # initial range
    alpha_q = min(q_aa, q_ab, q_ba, q_bb)
    beta_q = max(q_aa, q_ab, q_ba, q_bb)

    # analyze cubic boundary polynomials q1 to q4
    for i in 1:4
        c0 = 0.0
        c1 = 0.0
        c2 = 0.0
        c3 = 0.0

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
        else # i == 4
            c0 = c00 + r*c10 + r^2*c20 + r^3*c30
            c1 = c01 + r*c11 + r^2*c21
            c2 = c02 + r*c12
            c3 = c03
        end

        if c3 == 0
            if abs(c1) < 2*abs(c2)*r
                e = c0 - c1^2/(4*c2)
                if c2 > 0
                    alpha_q = min(alpha_q, e)
                else
                    beta_q = max(beta_q, e)
                end
            end
        else 
            delta = c2^2 - 3*c1*c3
            if delta > 0
                L = sign(c3) * (c1 + 3*c3*r^2)
                R = 2*abs(c2)*r
                if L > R
                    if abs(c2) < 3*abs(c3)*r
                        xm = -(c2 - sqrt(delta)) / (3*c3)
                        xp = -(c2 + sqrt(delta)) / (3*c3)
                        alpha_q = min(alpha_q, ((c3*xm + c2)*xm + c1)*xm + c0)
                        beta_q  = max(beta_q, ((c3*xp + c2)*xp + c1)*xp + c0)
                    end
                elseif L > -R 
                    if c2 > 0
                        xm = -(c2 - sqrt(delta)) / (3*c3)
                        alpha_q = min(alpha_q, ((c3*xm + c2)*xm + c1)*xm + c0)
                    elseif c2 < 0
                        xp = -(c2 + sqrt(delta)) / (3*c3)
                        beta_q = max(beta_q, ((c3*xp + c2)*xp + c1)*xp + c0)
                    end
                end
            end
        end
    end

    alpha_q, beta_q = interior!(
        alpha_q, beta_q,
        c00,c10,c01,c20,c11,c02,c30,c21,c12,c03,
        r,
    )

    r_aa =  r^4*(c31 + c22 + c13) - r^5*(c32 + c23) + r^6*c33
    r_ab = -r^4*(c31 - c22 + c13) - r^5*(c32 - c23) - r^6*c33
    r_ba = -r^4*(c31 - c22 + c13) + r^5*(c32 - c23) - r^6*c33
    r_bb =  r^4*(c31 + c22 + c13) + r^5*(c32 + c23) + r^6*c33

    alpha_r = minimum((r_aa, r_ab, r_ba, r_bb))
    beta_r  = maximum((r_aa, r_ab, r_ba, r_bb))

    for i in 1:4
        c0 = 0.0
        c1 = 0.0
        c2 = 0.0
        c3 = 0.0

        if i == 1
            c0 = 0
            c1 = -c13*r^3
            c2 =  c22*r^2 - c23*r^3
            c3 = -c31*r + c32*r^2 - c33*r^3
        elseif i == 2
            c0 = 0
            c1 =  c13*r^3
            c2 =  c22*r^2 + c23*r^3
            c3 =  c31*r + c32*r^2 + c33*r^3
        elseif i == 3
            c0 = 0
            c1 = -c31*r^3
            c2 =  c22*r^2 - c32*r^3
            c3 = -c13*r + c23*r^2 - c33*r^3
        else # i == 4
            c0 = 0
            c1 =  c31*r^3
            c2 =  c22*r^2 + c32*r^3
            c3 =  c13*r + c23*r^2 + c33*r^3
        end

        if c3 == 0
            if abs(c1) < 2*abs(c2)*r
                e = c0 - c1^2 / (4*c2)
                if c2 > 0
                    alpha_r = min(alpha_r, e)
                else
                    beta_r  = max(beta_r, e)
                end
            end
        else
            delta = c2^2 - 3*c1*c3
            if delta > 0
                L = sign(c3) * (c1 + 3*c3*r^2)
                R = 2*abs(c2)*r
                if L > R
                    if abs(c2) < 3*abs(c3)*r
                        xm = -(c2 - sqrt(delta)) / (3*c3)
                        xp = -(c2 + sqrt(delta)) / (3*c3)
                        alpha_r = min(alpha_r, ((c3*xm + c2)*xm + c1)*xm + c0)
                        beta_r  = max(beta_r, ((c3*xp + c2)*xp + c1)*xp + c0)
                    end
                elseif L > -R
                    if c2 > 0
                        xm = -(c2 - sqrt(delta)) / (3*c3)
                        alpha_r = min(alpha_r, ((c3*xm + c2)*xm + c1)*xm + c0)
                    elseif c2 < 0
                        xp = -(c2 + sqrt(delta)) / (3*c3)
                        beta_r  = max(beta_r, ((c3*xp + c2)*xp + c1)*xp + c0)
                    end
                end
            end
        end
    end

    alpha_r, beta_r = update_alpha_beta_r_interior!(alpha_r, beta_r,
    c31,c22,c13,c32,c23,c33, r)

    return alpha_q + alpha_r, beta_q + beta_r
end




#=
    hermite4(f, df, B, d)

Range of `f` over box `B = (mx, my, r)` approximated by the
Hermite form of order 4 (Section 4 and Appendix B.5 in the paper).

Arguments:
- poly_str::String: the original scalar function as a string
- B = (mx, my, r): center (mx,my) and half-width r of the box
- d::Integer:     total degree of f(x, y)

Returns:
    (lower_bound, upper_bound)
=#
function hermite4(poly::Polynomial, B::myBox, degree::Int; sharing::Bool=true)::myInterval
    global max_x, max_y

    mx, my = B.x_2, B.y_2
    rx, ry = B.rx, B.ry

    r = rx

    n = fld(degree, 4)   # floor(d/4)

    trinomial_table!(n)

    evaluate_hermite(poly, B, sharing)

    # -------------
    # 1. Hermite bicubic for f
    # -------------
    # Arrays for F, Fx, Fy, Fxy at the four corners
    # Maple used indices 0,1, so we shift to 1,2.
    F   = zeros(Float64, 2, 2)
    Fx  = zeros(Float64, 2, 2)
    Fy  = zeros(Float64, 2, 2)
    Fxy = zeros(Float64, 2, 2)

    for i in 0:1, j in 0:1
        ix = (i == 0) ? 1 : 3
        iy = (j == 0) ? 1 : 3
        pt = get_point(B, ix, iy)

        F[i+1,  j+1]  = derivatives_hermite[pt][1][1]
        Fx[i+1, j+1]  = derivatives_hermite[pt][2][1]
        Fy[i+1, j+1]  = derivatives_hermite[pt][1][2]
        Fxy[i+1, j+1] = derivatives_hermite[pt][2][2]
    end

    # Coefficients of Hermite bicubic p in centered form
    c00 = ( Fxy[1,1] - Fxy[1,2] - Fxy[2,1] + Fxy[2,2])*r^2/16 +
          (  Fx[1,1] + Fx[1,2] - Fx[2,1] - Fx[2,2])*r/8   +
          (  Fy[1,1] - Fy[1,2] + Fy[2,1] - Fy[2,2])*r/8   +
          (   F[1,1] + F[1,2] + F[2,1] + F[2,2])/4

    c10 = (-Fxy[1,1] + Fxy[1,2] - Fxy[2,1] + Fxy[2,2])*r/16 +
          (- Fx[1,1] - Fx[1,2] - Fx[2,1] - Fx[2,2])/8      +
          (- Fy[1,1] + Fy[1,2] + Fy[2,1] - Fy[2,2])*3/16   +
          (-  F[1,1] - F[1,2] + F[2,1] + F[2,2])*3/(8r)

    c01 = (-Fxy[1,1] - Fxy[1,2] + Fxy[2,1] + Fxy[2,2])*r/16 +
          (- Fx[1,1] + Fx[1,2] + Fx[2,1] - Fx[2,2])*3/16    +
          (- Fy[1,1] - Fy[1,2] - Fy[2,1] - Fy[2,2])/8       +
          (-  F[1,1] + F[1,2] - F[2,1] + F[2,2])*3/(8r)

    c20 = (-Fxy[1,1] + Fxy[1,2] + Fxy[2,1] - Fxy[2,2])/16 +
          (- Fx[1,1] - Fx[1,2] + Fx[2,1] + Fx[2,2])/(8r)

    c11 = ( Fxy[1,1] + Fxy[1,2] + Fxy[2,1] + Fxy[2,2])/16 +
          (  Fx[1,1] - Fx[1,2] + Fx[2,1] - Fx[2,2])*3/(16r) +
          (  Fy[1,1] + Fy[1,2] - Fy[2,1] - Fy[2,2])*3/(16r) +
          (   F[1,1] - F[1,2] - F[2,1] + F[2,2])*9/(16r^2)

    c02 = (-Fxy[1,1] + Fxy[1,2] + Fxy[2,1] - Fxy[2,2])/16 +
          (- Fy[1,1] + Fy[1,2] - Fy[2,1] + Fy[2,2])/(8r)

    c30 = ( Fxy[1,1] - Fxy[1,2] + Fxy[2,1] - Fxy[2,2])/(16r) +
          (  Fx[1,1] + Fx[1,2] + Fx[2,1] + Fx[2,2])/(8r^2)  +
          (  Fy[1,1] - Fy[1,2] - Fy[2,1] + Fy[2,2])/(16r^2) +
          (   F[1,1] + F[1,2] - F[2,1] - F[2,2])/(8r^3)

    c21 = ( Fxy[1,1] + Fxy[1,2] - Fxy[2,1] - Fxy[2,2])/(16r) +
          (  Fx[1,1] - Fx[1,2] - Fx[2,1] + Fx[2,2])*3/(16r^2)

    c12 = ( Fxy[1,1] - Fxy[1,2] + Fxy[2,1] - Fxy[2,2])/(16r) +
          (  Fy[1,1] - Fy[1,2] - Fy[2,1] + Fy[2,2])*3/(16r^2)

    c03 = ( Fxy[1,1] + Fxy[1,2] - Fxy[2,1] - Fxy[2,2])/(16r) +
          (  Fx[1,1] - Fx[1,2] - Fx[2,1] + Fx[2,2])/(16r^2) +
          (  Fy[1,1] + Fy[1,2] + Fy[2,1] + Fy[2,2])/(8r^2)   +
          (   F[1,1] - F[1,2] + F[2,1] - F[2,2])/(8r^3)

    c31 = (-Fxy[1,1] - Fxy[1,2] - Fxy[2,1] - Fxy[2,2])/(16r^2) +
          (- Fx[1,1] + Fx[1,2] - Fx[2,1] + Fx[2,2])*3/(16r^3) +
          (- Fy[1,1] - Fy[1,2] + Fy[2,1] + Fy[2,2])/(16r^3)   +
          (-  F[1,1] + F[1,2] + F[2,1] - F[2,2])*3/(16r^4)

    c22 = ( Fxy[1,1] - Fxy[1,2] - Fxy[2,1] + Fxy[2,2])/(16r^2)

    c13 = (-Fxy[1,1] - Fxy[1,2] - Fxy[2,1] - Fxy[2,2])/(16r^2) +
          (- Fx[1,1] + Fx[1,2] - Fx[2,1] + Fx[2,2])/(16r^3)   +
          (- Fy[1,1] - Fy[1,2] + Fy[2,1] + Fy[2,2])*3/(16r^3) +
          (-  F[1,1] + F[1,2] + F[2,1] - F[2,2])*3/(16r^4)

    c32 = (-Fxy[1,1] + Fxy[1,2] - Fxy[2,1] + Fxy[2,2])/(16r^3) +
          (- Fy[1,1] + Fy[1,2] + Fy[2,1] - Fy[2,2])/(16r^4)

    c23 = (-Fxy[1,1] - Fxy[1,2] + Fxy[2,1] + Fxy[2,2])/(16r^3) +
          (- Fx[1,1] + Fx[1,2] + Fx[2,1] - Fx[2,2])/(16r^4)

    c33 = ( Fxy[1,1] + Fxy[1,2] + Fxy[2,1] + Fxy[2,2])/(16r^4) +
          (  Fx[1,1] - Fx[1,2] + Fx[2,1] - Fx[2,2])/(16r^5)   +
          (  Fy[1,1] + Fy[1,2] - Fy[2,1] - Fy[2,2])/(16r^5)   +
          (   F[1,1] - F[1,2] - F[2,1] + F[2,2])/(16r^6)

    # range of p using split form
    a0, b0 = range_bi3(
        c00,c10,c01,c20,c11,c02,c30,c21,c12,c03,
        c31,c22,c13,c32,c23,c33,
        r,
    )

    # -------------
    # 2. Add remainder bound
    # -------------
    omega = (1/24) * r^4
    V4 = 0.0

    for k in 1:n
        for j in 0:k
            ub = k - j
            vb = j

            xi0 = 4*ub 
            xi1 = 4*ub + 1
            yi0 = 4*vb
            yi1 = 4*vb + 1

            for ii in 0:1, jj in 0:1
                ix = (ii == 0) ? 1 : 3
                iy = (jj == 0) ? 1 : 3
                pt  = get_point(B, ix, iy)

                tab = derivatives_hermite[pt]

                F[ii+1,  jj+1] = get_cached_deriv(tab, xi0, yi0, max_x, max_y, degree)
                Fx[ii+1, jj+1] = get_cached_deriv(tab, xi1, yi0, max_x, max_y, degree)
                Fy[ii+1, jj+1] = get_cached_deriv(tab, xi0, yi1, max_x, max_y, degree)
                Fxy[ii+1,jj+1] = get_cached_deriv(tab, xi1, yi1, max_x, max_y, degree)

                # F[ii+1,  jj+1]  = derivatives_hermite[pt][xi0][yi0] 
                # Fx[ii+1, jj+1]  = derivatives_hermite[pt][xi1][yi0]  
                # Fy[ii+1, jj+1]  = derivatives_hermite[pt][xi0][yi1]  
                # Fxy[ii+1,jj+1]  = derivatives_hermite[pt][xi1][yi1]  
            end

            # same coefficient formulas, now applied to these derivatives:
            c00 = ( Fxy[1,1] - Fxy[1,2] - Fxy[2,1] + Fxy[2,2])*r^2/16 +
                  (  Fx[1,1] + Fx[1,2] - Fx[2,1] - Fx[2,2])*r/8   +
                  (  Fy[1,1] - Fy[1,2] + Fy[2,1] - Fy[2,2])*r/8   +
                  (   F[1,1] + F[1,2] + F[2,1] + F[2,2])/4

            c10 = (-Fxy[1,1] + Fxy[1,2] - Fxy[2,1] + Fxy[2,2])*r/16 +
                  (- Fx[1,1] - Fx[1,2] - Fx[2,1] - Fx[2,2])/8      +
                  (- Fy[1,1] + Fy[1,2] + Fy[2,1] - Fy[2,2])*3/16   +
                  (-  F[1,1] - F[1,2] + F[2,1] + F[2,2])*3/(8r)

            c01 = (-Fxy[1,1] - Fxy[1,2] + Fxy[2,1] + Fxy[2,2])*r/16 +
                  (- Fx[1,1] + Fx[1,2] + Fx[2,1] - Fx[2,2])*3/16    +
                  (- Fy[1,1] - Fy[1,2] - Fy[2,1] - Fy[2,2])/8       +
                  (-  F[1,1] + F[1,2] - F[2,1] + F[2,2])*3/(8r)

            c20 = (-Fxy[1,1] + Fxy[1,2] + Fxy[2,1] - Fxy[2,2])/16 +
                  (- Fx[1,1] - Fx[1,2] + Fx[2,1] + Fx[2,2])/(8r)

            c11 = ( Fxy[1,1] + Fxy[1,2] + Fxy[2,1] + Fxy[2,2])/16 +
                  (  Fx[1,1] - Fx[1,2] + Fx[2,1] - Fx[2,2])*3/(16r) +
                  (  Fy[1,1] + Fy[1,2] - Fy[2,1] - Fy[2,2])*3/(16r) +
                  (   F[1,1] - F[1,2] - F[2,1] + F[2,2])*9/(16r^2)

            c02 = (-Fxy[1,1] + Fxy[1,2] + Fxy[2,1] - Fxy[2,2])/16 +
                  (- Fy[1,1] + Fy[1,2] - Fy[2,1] + Fy[2,2])/(8r)

            c30 = ( Fxy[1,1] - Fxy[1,2] + Fxy[2,1] - Fxy[2,2])/(16r) +
                  (  Fx[1,1] + Fx[1,2] + Fx[2,1] + Fx[2,2])/(8r^2)  +
                  (  Fy[1,1] - Fy[1,2] - Fy[2,1] + Fy[2,2])/(16r^2) +
                  (   F[1,1] + F[1,2] - F[2,1] - F[2,2])/(8r^3)

            c21 = ( Fxy[1,1] + Fxy[1,2] - Fxy[2,1] - Fxy[2,2])/(16r) +
                  (  Fx[1,1] - Fx[1,2] - Fx[2,1] + Fx[2,2])*3/(16r^2)

            c12 = ( Fxy[1,1] - Fxy[1,2] + Fxy[2,1] - Fxy[2,2])/(16r) +
                  (  Fy[1,1] - Fy[1,2] - Fy[2,1] + Fy[2,2])*3/(16r^2)

            c03 = ( Fxy[1,1] + Fxy[1,2] - Fxy[2,1] - Fxy[2,2])/(16r) +
                  (  Fx[1,1] - Fx[1,2] - Fx[2,1] + Fx[2,2])/(16r^2) +
                  (  Fy[1,1] + Fy[1,2] + Fy[2,1] + Fy[2,2])/(8r^2)   +
                  (   F[1,1] - F[1,2] + F[2,1] - F[2,2])/(8r^3)

            c31 = (-Fxy[1,1] - Fxy[1,2] - Fxy[2,1] - Fxy[2,2])/(16r^2) +
                  (- Fx[1,1] + Fx[1,2] - Fx[2,1] + Fx[2,2])*3/(16r^3) +
                  (- Fy[1,1] - Fy[1,2] + Fy[2,1] + Fy[2,2])/(16r^3)   +
                  (-  F[1,1] + F[1,2] + F[2,1] - F[2,2])*3/(16r^4)

            c22 = ( Fxy[1,1] - Fxy[1,2] - Fxy[2,1] + Fxy[2,2])/(16r^2)

            c13 = (-Fxy[1,1] - Fxy[1,2] - Fxy[2,1] - Fxy[2,2])/(16r^2) +
                  (- Fx[1,1] + Fx[1,2] - Fx[2,1] + Fx[2,2])/(16r^3)   +
                  (- Fy[1,1] - Fy[1,2] + Fy[2,1] + Fy[2,2])*3/(16r^3) +
                  (-  F[1,1] + F[1,2] + F[2,1] - F[2,2])*3/(16r^4)

            c32 = (-Fxy[1,1] + Fxy[1,2] - Fxy[2,1] + Fxy[2,2])/(16r^3) +
                  (- Fy[1,1] + Fy[1,2] + Fy[2,1] - Fy[2,2])/(16r^4)

            c23 = (-Fxy[1,1] - Fxy[1,2] + Fxy[2,1] + Fxy[2,2])/(16r^3) +
                  (- Fx[1,1] + Fx[1,2] + Fx[2,1] - Fx[2,2])/(16r^4)

            c33 = ( Fxy[1,1] + Fxy[1,2] + Fxy[2,1] + Fxy[2,2])/(16r^4) +
                  (  Fx[1,1] - Fx[1,2] + Fx[2,1] - Fx[2,2])/(16r^5)   +
                  (  Fy[1,1] + Fy[1,2] - Fy[2,1] - Fy[2,2])/(16r^5)   +
                  (   F[1,1] - F[1,2] - F[2,1] + F[2,2])/(16r^6)

            va, vb = range_bi3(
                c00,c10,c01,c20,c11,c02,c30,c21,c12,c03,
                c31,c22,c13,c32,c23,c33,
                r,
            )

            V4 += omega^k * trinomial(k, j) * max(abs(va), abs(vb))
        end
    end

    # Maple had simplify; we just return the numeric quantities
    return myInterval(a0 - V4, b0 + V4)
end

# B = (0.0, 0.0, 1.0)
# poly_str = "x^4 + 2x^2*y^2 + y^4"
# d = 4
# lo, hi = hermite4(poly_str, B, d)
# println("Range ≈ [$lo, $hi]")
