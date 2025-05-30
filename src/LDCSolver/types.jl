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

struct ConvectionOperator{T} <: AbstractNonLinearOperator where T<:Number
    Muu::AbstractMatrix{T}
    Muv::AbstractMatrix{T}
    Mvu::AbstractMatrix{T}
    Mvv::AbstractMatrix{T}
    function ConvectionOperator(mesh::CartesianMesh, ::Type{T}=Float64) where T<:Number
        Muu = spzeros(T, mesh.N^2, mesh.N^2)
        Muv = spzeros(T, mesh.N^2, mesh.N^2)
        Mvu = spzeros(T, mesh.N^2, mesh.N^2)
        Mvv = spzeros(T, mesh.N^2, mesh.N^2)
        return new{T}(Muu, Muv, Mvu, Mvv)
    end
end
function ConvectionOperator(mesh::CartesianMesh)
    return ConvectionOperator(mesh, Float64)
end

function update_convection_operators!(op::ConvectionOperator,
    mesh::CartesianMesh,
    uvp,
    limiter::AbstractLimiter,
    bc_u::DomainBoundaryConditions,
    bc_v::DomainBoundaryConditions,
)
    any(isnan.(uvp)) && @warn "uvp contains nans"
    N = mesh.N
    @assert length(uvp) == 3(N^2)
    u = @view uvp[1:N^2]
    v = @view uvp[N^2+1:2N^2]
    set_convection_operators!(op.Muu, op.Muv, mesh, u, u, v, limiter, bc_u)
    set_convection_operators!(op.Mvu, op.Mvv, mesh, v, u, v, limiter, bc_v)
    nothing
end

# Helper function to create ConvectionOperator with element type matching uvp
function create_convection_operator(mesh::CartesianMesh, uvp::AbstractVector{T}) where T
    return ConvectionOperator(mesh, T)
end

struct LDCOperators{T<:Number}
    convection::ConvectionOperator{T}
    dx_pressure::LinearOperator
    dy_pressure::LinearOperator
    laplacian_u::LinearOperator
    laplacian_v::LinearOperator
    continuity_u::LinearOperator
    continuity_v::LinearOperator
    function LDCOperators(mesh::CartesianMesh, u_bc::DomainBoundaryConditions, v_bc::DomainBoundaryConditions, ::Type{T}=Float64) where T<:Number
        convection = ConvectionOperator(mesh, T)
        dx_pressure = LinearOperator(get_pressure_dx_operator(mesh)...)
        dy_pressure = LinearOperator(get_pressure_dy_operator(mesh)...)
        laplacian_u = LinearOperator(get_laplacian_operator(mesh, u_bc.north.value)...)
        laplacian_v = LinearOperator(get_laplacian_operator(mesh, v_bc.north.value)...)

        continuity_u = LinearOperator(get_continuity_operator_x(mesh)...)
        continuity_v = LinearOperator(get_continuity_operator_y(mesh)...)

        return new{T}(convection, dx_pressure, dy_pressure, laplacian_u, laplacian_v, continuity_u, continuity_v)
    end
end

abstract type AbstractLDCProblem end

struct LDCProblem{T<:Number} <: AbstractLDCProblem
    mesh::CartesianMesh
    ν::Float64
    ψ::AbstractLimiter
    bc_u::DomainBoundaryConditions
    bc_v::DomainBoundaryConditions
    ops::LDCOperators{T}
    function LDCProblem(mesh::CartesianMesh, ν::Float64, ψ::AbstractLimiter, U=1.0, ::Type{T}=Float64) where T<:Number
        bc_u = DomainBoundaryConditions(north=DirichletBC(U), south=DirichletBC(0.0), east=DirichletBC(0.0), west=DirichletBC(0.0))
        bc_v = DomainBoundaryConditions(north=DirichletBC(0.0), south=DirichletBC(0.0), east=DirichletBC(0.0), west=DirichletBC(0.0))

        ops = LDCOperators(mesh, bc_u, bc_v, T)
        return new{T}(mesh, ν, ψ, bc_u, bc_v, ops)
    end
end

