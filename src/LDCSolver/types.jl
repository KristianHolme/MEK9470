
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
    Muu::AbstractMatrix
    Muv::AbstractMatrix
    Mvu::AbstractMatrix
    Mvv::AbstractMatrix
    function ConvectionOperator(mesh::CartesianMesh)
        Muu = spzeros(Float64, mesh.N^2, mesh.N^2)
        Muv = spzeros(Float64, mesh.N^2, mesh.N^2)
        Mvu = spzeros(Float64, mesh.N^2, mesh.N^2)
        Mvv = spzeros(Float64, mesh.N^2, mesh.N^2)
        return new(Muu, Muv, Mvu, Mvv)
    end
end

function update_convection_operators!(op::ConvectionOperator,
    mesh::CartesianMesh,
    uvp,
    limiter::AbstractLimiter,
    bc_u::DomainBoundaryConditions,
    bc_v::DomainBoundaryConditions,
)
    N = mesh.N
    @assert length(uvp) == 3(N^2)
    u = @view uvp[1:N^2]
    v = @view uvp[N^2+1:2N^2]
    set_convection_operators!(op.Muu, op.Muv, mesh, u, u, v, limiter, bc_u)
    set_convection_operators!(op.Mvu, op.Mvv, mesh, v, u, v, limiter, bc_v)
    nothing
end


struct LDCOperators
    convection::ConvectionOperator
    dx_pressure::LinearOperator
    dy_pressure::LinearOperator
    laplacian_u::LinearOperator
    laplacian_v::LinearOperator
    continuity_u::LinearOperator
    continuity_v::LinearOperator
    function LDCOperators(mesh::CartesianMesh, u_bc::DomainBoundaryConditions, v_bc::DomainBoundaryConditions)
        convection = ConvectionOperator(mesh)
        dx_pressure = LinearOperator(get_pressure_dx_operator(mesh)...)
        dy_pressure = LinearOperator(get_pressure_dy_operator(mesh)...)
        laplacian_u = LinearOperator(get_laplacian_operator(mesh, u_bc.north.value)...)
        laplacian_v = LinearOperator(get_laplacian_operator(mesh, v_bc.north.value)...)

        continuity_u = LinearOperator(get_continuity_operator(mesh)...)

        Mcv, bcv = get_continuity_operator(mesh)
        continuity_v = LinearOperator(Mcv', bcv)
        return new(convection, dx_pressure, dy_pressure, laplacian_u, laplacian_v, continuity_u, continuity_v)
    end
end

abstract type AbstractLDCProblem end

struct LDCProblem <: AbstractLDCProblem
    mesh::CartesianMesh
    ν::Float64
    ψ::AbstractLimiter
    bc_u::DomainBoundaryConditions
    bc_v::DomainBoundaryConditions
    ops::LDCOperators
    function LDCProblem(mesh::CartesianMesh, ν::Float64, ψ::AbstractLimiter, U=1.0)
        bc_u = DomainBoundaryConditions(north=DirichletBC(U), south=DirichletBC(0.0), east=DirichletBC(0.0), west=DirichletBC(0.0))
        bc_v = DomainBoundaryConditions(north=DirichletBC(0.0), south=DirichletBC(0.0), east=DirichletBC(0.0), west=DirichletBC(0.0))

        ops = LDCOperators(mesh, bc_u, bc_v)
        return new(mesh, ν, ψ, bc_u, bc_v, ops)
    end
end

