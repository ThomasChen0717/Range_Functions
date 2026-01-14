#= 
    @author: Thomas Chen
    @advisors: Yap Chee, Kai Hormann, Bingwei Zhang
    @date: 2025-06-03
    @description: A mutable struct representing a box region in a two-dimensional space.
=# 

#=
    myBox

    A mutable struct representing a box region in a two-dimensional space.

    # Fields
    - `x`: The x interval of the box
    - `y`: The y interval of the box
    - `depth`: The depth of the box
    - `up`: A vector of boxes representing the upper neighbors
    - `down`: A vector of boxes representing the lower neighbors
    - `left`: A vector of boxes representing the left neighbors
    - `right`: A vector of boxes representing the right neighbors
    - `x_1`: The lower bound of the x interval
    - `x_2`: The midpoint of the x interval
    - `x_3`: The upper bound of the x interval
    - `rx`: The radius of the x interval
    - `y_1`: The lower bound of the y interval
    - `y_2`: The midpoint of the y interval
    - `y_3`: The upper bound of the y interval
    - `ry`: The radius of the y interval
    - `pts_matrix`: A matrix containing the coordinates of the box points
    - `FB`: An array of intervals representing the function values at the box points. Nothing represents uninitialized 
    - `QB`: An array of vectors of functions representing the derivative values at the box points
=#
mutable struct myBox
    x::myInterval
    y::myInterval
    depth::Int 
    up::Vector{myBox}
    down::Vector{myBox}
    left::Vector{myBox}
    right::Vector{myBox}
    x_1::Float64
    x_2::Float64
    x_3::Float64
    rx::Float64
    y_1::Float64
    y_2::Float64
    y_3::Float64
    ry::Float64
    pts_matrix::Array{Tuple{Float64, Float64}}
    FB::Array{Union{myInterval, Nothing}}
    QB::Array{Vector{Union{Float64, Function}}}
end

#=
    myBox(x::myInterval, y::myInterval, poly::Polynomial; depth::Int=0)::myBox

    Constructor for the `myBox` struct. Initializes a new box with the given x and y intervals, depth, and default values for other attributes.

    # Arguments
    - `x`: The x interval of the box
    - `y`: The y interval of the box
    - `poly`: The polynomial of the box
    - `depth`: The depth of the box (default: 0)

    # Returns
    - `myBox`: An instance of the `myBox` struct
=#
function myBox(x::myInterval, y::myInterval; depth::Int=0)::myBox
    x_1 = x.lower
    x_3 = x.upper
    x_2 = (x_1 + x_3) / 2
    rx = (x_3 - x_1) / 2
    
    y_1 = y.lower
    y_3 = y.upper
    y_2 = (y_1 + y_3) / 2
    ry = (y_3 - y_1) / 2

    pts_matrix = Array{Tuple{Float64, Float64}}(undef, 3, 3)

    pts_matrix[1, 1] = (x_1, y_1)  
    pts_matrix[1, 2] = (x_1, y_2)
    pts_matrix[1, 3] = (x_1, y_3)
    
    pts_matrix[2, 1] = (x_2, y_1) 
    pts_matrix[2, 2] = (x_2, y_2) 
    pts_matrix[2, 3] = (x_2, y_3)  
    
    pts_matrix[3, 1] = (x_3, y_1)  
    pts_matrix[3, 2] = (x_3, y_2)  
    pts_matrix[3, 3] = (x_3, y_3) 

    max_x_multiples = div(max_x, 3) + 1
    max_y_multiples = div(max_y, 3) + 1

    FB = Array{Union{myInterval, Nothing}}(undef, max_x_multiples + 1, max_y_multiples + 1)
    fill!(FB, nothing)

    QB = Array{Vector{Union{Float64, Function}}}(undef, max_x_multiples + 1, max_y_multiples + 1)
    for i in 1:(max_x_multiples + 1), j in 1:(max_y_multiples + 1)
        QB[i,j] = []
    end

    return myBox(x, y, depth, [],[],[],[], x_1, x_2, x_3,rx, y_1, y_2, y_3, ry, pts_matrix, FB, QB)
end


#=
    internal_neighbors!(children::Vector{myBox})

    Function to establish internal neighbors between the given children boxes.

    # Arguments
    - `children`: A vector of four boxes representing the children boxes
=#
function internal_neighbors!(children::Vector{myBox})
    # children order: [bottom-left, top-left, bottom-right, top-right]
    bottom_left, top_left, bottom_right, top_right = children

    push!(bottom_left.up, top_left)     
    push!(bottom_left.right, bottom_right) 

    
    push!(top_left.down, bottom_left)    
    push!(top_left.right, top_right)     

    
    push!(bottom_right.up, top_right)    
    push!(bottom_right.left, bottom_left) 

    
    push!(top_right.down, bottom_right)  
    push!(top_right.left, top_left)      
end

#=
    external_neighbors!(parent::myBox, children::Vector{myBox})

    Function to establish external neighbors between the given parent box and its children boxes.

    # Arguments
    - `parent`: The parent box
    - `children`: A vector of four boxes representing the children boxes
=#
function external_neighbors!(parent::myBox, children::Vector{myBox})
    bottom_left, top_left, bottom_right, top_right = children

    for up_neighbor in parent.up
        filter!(box -> box !== parent, up_neighbor.down)

        if x_overlap(up_neighbor, top_left)
            push!(up_neighbor.down, top_left)
            push!(top_left.up, up_neighbor)
        end

        if x_overlap(up_neighbor, top_right)
            push!(up_neighbor.down, top_right)
            push!(top_right.up, up_neighbor)
        end
    end

    for down_neighbor in parent.down
        filter!(box -> box!== parent, down_neighbor.up)
        if x_overlap(down_neighbor, bottom_left)
            push!(down_neighbor.up, bottom_left)
            push!(bottom_left.down, down_neighbor)
        end
        if x_overlap(down_neighbor, bottom_right)
            push!(down_neighbor.up, bottom_right)
            push!(bottom_right.down, down_neighbor)
        end
    end
    for left_neighbor in parent.left
        filter!(box -> box!== parent, left_neighbor.right)
        if y_overlap(left_neighbor, bottom_left)
            push!(left_neighbor.right, bottom_left)
            push!(bottom_left.left, left_neighbor)
        end
        if y_overlap(left_neighbor, top_left)
            push!(left_neighbor.right, top_left)
            push!(top_left.left, left_neighbor)
        end
    end
    for right_neighbor in parent.right
        filter!(box -> box!== parent, right_neighbor.left)
        if y_overlap(right_neighbor, bottom_right)
            push!(right_neighbor.left, bottom_right)
            push!(bottom_right.right, right_neighbor)
        end
        if y_overlap(right_neighbor, top_right)
            push!(right_neighbor.left, top_right)
            push!(top_right.right, right_neighbor)
        end
    end
end

#= 
    split_box(b::myBox, polynomial::Polynomial)::Vector{myBox}

    Function to split the given box into four equal sub-boxes.

    # Arguments
    - `b`: The box to split
    - `polynomial`: The polynomial of the box

    # Returns
    - `Vector{myBox}`: A vector of four sub-boxes
=#
function split_box(b::myBox)::Vector{myBox}
    mx = (b.x.lower + b.x.upper) / 2
    my = (b.y.lower + b.y.upper) / 2

    x1 = myInterval(b.x.lower, mx)
    x2 = myInterval(mx, b.x.upper)
    y1 = myInterval(b.y.lower, my)
    y2 = myInterval(my, b.y.upper)

    d  = b.depth + 1

    children = [
        myBox(x1, y1; depth=d),
        myBox(x1, y2; depth=d),
        myBox(x2, y1; depth=d),
        myBox(x2, y2; depth=d)
    ]

    internal_neighbors!(children)

    external_neighbors!(b, children)

    return children
end



#=
    y_overlap(b::myBox, c::myBox)::Bool

    Function to check if the y-intervals of two boxes overlap.

    # Arguments
    - `b`: The first box
    - `c`: The second box

    # Returns
    - `Bool`: `true` if the y-intervals overlap, `false` otherwise
=#
function y_overlap(b::myBox, c::myBox)::Bool
    return (b.y.lower < c.y.upper) && (b.y.upper > c.y.lower)
end

#=
    x_overlap(b::myBox, c::myBox)::Bool

    Function to check if the x-intervals of two boxes overlap.

    # Arguments
    - `b`: The first box
    - `c`: The second box

    # Returns
    - `Bool`: `true` if the x-intervals overlap, `false` otherwise
=#
function x_overlap(b::myBox, c::myBox)::Bool
    return (b.x.lower < c.x.upper) && (b.x.upper > c.x.lower)
end







