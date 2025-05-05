
abstract type AbstractDiscretization end
abstract type AbstractLimiter end
abstract type AbstractMesh end

struct CartesianMesh <: AbstractMesh
    N::Int
end

abstract type AbstractLDCProblem end

struct LDCProblem <: AbstractLDCProblem
    mesh::CartesianMesh
    f::Function
    ν::Float64
    limiter::AbstractLimiter
end


"""
LinearOperator

a linear numerical operator approximation of the continuous operator f, such that
f(u) ≈ M*u + b
"""

abstract type AbstractOperator end
abstract type AbstractNonLinearOperator <: AbstractOperator end
struct LinearOperator <: AbstractOperator
    M::AbstractMatrix
    b::AbstractVector
end

struct ConvectionOperator <: AbstractNonLinearOperator
    M::AbstractMatrix
    b::AbstractVector
end

function update_matrix! end

struct LDCOperators
    convection::ConvectionOperator
    dx_pressure::LinearOperator
    dy_pressure::LinearOperator
    laplacian_u::LinearOperator
    laplacian_v::LinearOperator
    continuity::LinearOperator
end