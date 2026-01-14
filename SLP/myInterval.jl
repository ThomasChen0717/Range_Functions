#= 
    @author: Thomas Chen
    @advisors: Yap Chee, Kai Hormann, Bingwei Zhang
    @date: 2025-06-17
    @description: myInterval Struct and Operations
=# 

# interval.jl
struct myInterval
    lower::Float64
    upper::Float64
end

# myInterval addition
Base.:+(a::myInterval, b::myInterval) = myInterval(a.lower + b.lower, a.upper + b.upper)
Base.:+(a::Real, b::myInterval) = myInterval(a + b.lower, a + b.upper)
Base.:+(a::myInterval, b::Real) = myInterval(a.lower + b, a.upper + b)

# myInterval subtraction
Base.:-(a::myInterval, b::myInterval) = myInterval(a.lower - b.upper, a.upper - b.lower)
Base.:-(a::Real, b::myInterval) = myInterval(a - b.upper, a - b.lower)  # Fixed: swap b.upper and b.lower
Base.:-(a::myInterval, b::Real) = myInterval(a.lower - b, a.upper - b)


# myInterval multiplication
# More efficient interval multiplication
Base.:*(a::myInterval, b::myInterval) = begin
    # Use more efficient bounds calculation
    if a.lower >= 0 && b.lower >= 0
        # Both positive
        myInterval(a.lower * b.lower, a.upper * b.upper)
    elseif a.upper <= 0 && b.upper <= 0
        # Both negative
        myInterval(a.upper * b.upper, a.lower * b.lower)
    else
        # General case - only compute when necessary
        products = [a.lower * b.lower, a.lower * b.upper, a.upper * b.lower, a.upper * b.upper]
        myInterval(minimum(products), maximum(products))
    end
end
Base.:*(a::Real, I::myInterval) = myInterval(a,a) * I
Base.:*(I::myInterval, a::Real) = myInterval(a,a) * I



# myInterval exponentiation (integer power)
function Base.:^(a::myInterval, n::Int)
    if n == 0
        return myInterval(1.0, 1.0)
    elseif n == 1
        return a
    elseif n > 0
        result = myInterval(1.0, 1.0)
        base = a
        exponent = n
        
        while exponent > 0
            if exponent & 1 == 1  
                result = result * base
            end
            base = base * base
            exponent >>= 1 
        end
        
        return result
    else
        error("Negative interval exponentiation not supported")
    end
end
