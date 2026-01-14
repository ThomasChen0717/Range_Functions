#= 
    @author: Thomas Chen
    @advisors: Yap Chee, Kai Hormann, Bingwei Zhang
    @date: 2025-06-17
    @description: AST Types 
=# 

# Abstract base type for all expression nodes in the Abstract Syntax Tree (AST)
# This serves as the parent type for all mathematical expression representations
abstract type Expr end

# Represents a variable node in the AST
# Contains a symbolic name (e.g., :x, :y, :z) that can be used in mathematical expressions
struct VarNode <: Expr
    name::Symbol  
end

# Represents a floating-point constant node in the AST
# Contains a numerical value that remains constant throughout expression evaluation
struct FloatNode <: Expr
    value::Float64  
end

# Represents an addition operation node in the AST
# Performs binary addition between two sub-expressions (left + right)
struct AddNode <: Expr
    left::Expr   
    right::Expr  
end

# Represents a subtraction operation node in the AST
# Performs binary subtraction between two sub-expressions (left - right)
struct MinusNode <: Expr
    left::Expr   
    right::Expr  
end

# Represents a multiplication operation node in the AST
# Performs binary multiplication between two sub-expressions (left * right)
struct MulNode <: Expr
    left::Expr   
    right::Expr 
end

# Represents a power operation node in the AST
# Raises a base expression to a constant floating-point exponent (base^exponent)
struct PowNode <: Expr
    base::Expr        
    exponent::Float64 
end



