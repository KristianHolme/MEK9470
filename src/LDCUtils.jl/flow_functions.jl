
"""
    compute_circulation(ψ::AbstractVector, mesh::CartesianMesh, path_cells::Vector{Tuple{Int,Int}})

Compute circulation around a closed path using the streamfunction.
Circulation Γ = ∮ V⃗ · dl⃗ = Ψ_final - Ψ_initial for a closed path.

# Arguments
- `ψ`: Streamfunction vector
- `mesh`: CartesianMesh
- `path_cells`: Vector of (i,j) cell coordinates defining the closed path

# Returns
- Circulation value (scalar)
"""
function compute_circulation(ψ::AbstractVector, mesh::CartesianMesh, path_cells::Vector{Tuple{Int,Int}})
    N = mesh.N
    Ψ = reshape(ψ, N, N)

    if isempty(path_cells)
        return 0.0
    end

    # For a closed path, circulation is the difference in streamfunction
    # between any two points (since ψ is constant along streamlines)
    start_i, start_j = path_cells[1]
    end_i, end_j = path_cells[end]

    return Ψ[end_i, end_j] - Ψ[start_i, start_j]
end

"""
    compute_vorticity(uvp::AbstractVector, mesh::CartesianMesh)

Compute vorticity ω = ∂v/∂x - ∂u/∂y using finite differences.

# Arguments
- `uvp`: Combined velocity-pressure vector [u; v; p]
- `mesh`: CartesianMesh defining the computational domain

# Returns
- `ω`: Vorticity as a vector in the same ordering as velocity components
"""
function compute_vorticity(uvp::AbstractVector, mesh::CartesianMesh)
    u, v, p = split_uvp(uvp, mesh)
    N = mesh.N
    dx = dy = 1 / N

    U = reshape(u, N, N)
    V = reshape(v, N, N)

    ω = zeros(eltype(uvp), N, N)

    # Compute ∂v/∂x - ∂u/∂y using central differences in interior
    for i in 2:N-1, j in 2:N-1
        ∂v_∂x = (V[i+1, j] - V[i-1, j]) / (2 * dx)
        ∂u_∂y = (U[i, j+1] - U[i, j-1]) / (2 * dy)
        ω[i, j] = ∂v_∂x - ∂u_∂y
    end

    # Handle boundaries with forward/backward differences
    for j in 2:N-1
        # Left boundary (i=1)
        ∂v_∂x = (-3 * V[1, j] + 4 * V[2, j] - V[3, j]) / (2 * dx)
        ∂u_∂y = (U[1, j+1] - U[1, j-1]) / (2 * dy)
        ω[1, j] = ∂v_∂x - ∂u_∂y

        # Right boundary (i=N)
        ∂v_∂x = (3 * V[N, j] - 4 * V[N-1, j] + V[N-2, j]) / (2 * dx)
        ∂u_∂y = (U[N, j+1] - U[N, j-1]) / (2 * dy)
        ω[N, j] = ∂v_∂x - ∂u_∂y
    end

    for i in 2:N-1
        # Bottom boundary (j=1)
        ∂v_∂x = (V[i+1, 1] - V[i-1, 1]) / (2 * dx)
        ∂u_∂y = (-3 * U[i, 1] + 4 * U[i, 2] - U[i, 3]) / (2 * dy)
        ω[i, 1] = ∂v_∂x - ∂u_∂y

        # Top boundary (j=N)
        ∂v_∂x = (V[i+1, N] - V[i-1, N]) / (2 * dx)
        ∂u_∂y = (3 * U[i, N] - 4 * U[i, N-1] + U[i, N-2]) / (2 * dy)
        ω[i, N] = ∂v_∂x - ∂u_∂y
    end

    # Corner points (use one-sided differences)
    # Bottom-left (1,1)
    ∂v_∂x = (-3 * V[1, 1] + 4 * V[2, 1] - V[3, 1]) / (2 * dx)
    ∂u_∂y = (-3 * U[1, 1] + 4 * U[1, 2] - U[1, 3]) / (2 * dy)
    ω[1, 1] = ∂v_∂x - ∂u_∂y

    # Bottom-right (N,1)
    ∂v_∂x = (3 * V[N, 1] - 4 * V[N-1, 1] + V[N-2, 1]) / (2 * dx)
    ∂u_∂y = (-3 * U[N, 1] + 4 * U[N, 2] - U[N, 3]) / (2 * dy)
    ω[N, 1] = ∂v_∂x - ∂u_∂y

    # Top-left (1,N)
    ∂v_∂x = (-3 * V[1, N] + 4 * V[2, N] - V[3, N]) / (2 * dx)
    ∂u_∂y = (3 * U[1, N] - 4 * U[1, N-1] + U[1, N-2]) / (2 * dy)
    ω[1, N] = ∂v_∂x - ∂u_∂y

    # Top-right (N,N)
    ∂v_∂x = (3 * V[N, N] - 4 * V[N-1, N] + V[N-2, N]) / (2 * dx)
    ∂u_∂y = (3 * U[N, N] - 4 * U[N, N-1] + U[N, N-2]) / (2 * dy)
    ω[N, N] = ∂v_∂x - ∂u_∂y

    return vec(ω)
end


function compute_streamfunction(uvp, mesh)
    u, v, p = split_uvp(uvp, mesh)
    N = mesh.N
    dx = dy = 1 / N
    ω = compute_vorticity(uvp, mesh)

    U = reshape(u, N, N)
    V = reshape(v, N, N)

    A, b = get_laplacian_operator(mesh, 0.0)

    ψ = A \ -ω
end


"""
    verify_flow_relations(uvp::AbstractVector, mesh::CartesianMesh; tolerance=1e-6)

Verify the mathematical relationships between velocity, streamfunction, and potential.

# Arguments
- `uvp`: Combined velocity-pressure vector [u; v; p]
- `mesh`: CartesianMesh defining the computational domain  
- `tolerance`: Tolerance for numerical verification

# Returns
- Named tuple with verification results and error metrics
"""
function verify_flow_relations(uvp::AbstractVector, mesh::CartesianMesh; tolerance::Real=1e-6)
    u, v, p = split_uvp(uvp, mesh)
    N = mesh.N
    dx = dy = 1 / N

    # Compute flow functions
    ψ = compute_streamfunction(uvp, mesh)

    Ψ = reshape(ψ, N, N)
    U = reshape(u, N, N)
    V = reshape(v, N, N)

    ψ_x_error = 0.0
    ψ_y_error = 0.0
    count = 0
    for i in 2:N-1, j in 2:N-1
        ∂ψ_∂x = (Ψ[i+1, j] - Ψ[i-1, j]) / (2 * dx)
        ∂ψ_∂y = (Ψ[i, j+1] - Ψ[i, j-1]) / (2 * dy)

        ψ_x_error += (∂ψ_∂x - V[i, j])^2
        ψ_y_error += (-∂ψ_∂y - U[i, j])^2

        count += 1
    end

    ψ_x_rms = sqrt(ψ_x_error / count)
    ψ_y_rms = sqrt(ψ_y_error / count)

    return (
        streamfunction_valid=(ψ_x_rms < tolerance && ψ_y_rms < tolerance),
        ψ_x_rms_error=ψ_x_rms,
        ψ_y_rms_error=ψ_y_rms,
        max_ψ_error=max(ψ_x_rms, ψ_y_rms)
    )
end




