function split_uvp(uvp, mesh::CartesianMesh)
    N = mesh.N
    @assert length(uvp) == 3 * N^2
    u = uvp[1:N^2]
    v = uvp[N^2+1:2*N^2]
    p = uvp[2*N^2+1:3*N^2]
    return u, v, p
end

function continuity_equation(u, v, mesh::CartesianMesh)
    N = mesh.N
    dx = 1 / N
    r = zeros(N^2)
    r .= continuity_x(u, mesh) + continuity_y(v, mesh)
    r .*= dx
    return r
end

"""
continuity in x direction(for u). Transpose for v in y-direction
"""
function get_continuity_operator(mesh::CartesianMesh, type::Type{T}=Float64) where T<:AbstractFloat
    N = mesh.N
    dx = 1 / N
    M = spzeros(T, N^2, N^2)
    b = spzeros(T, N^2)
    for i in 1:N, j in 1:N
        center = cartesian_to_linear_indices(i, j, mesh)
        east = i != N ? cartesian_to_linear_indices(i + 1, j, mesh) : nothing
        west = i != 1 ? cartesian_to_linear_indices(i - 1, j, mesh) : nothing

        if i == 1
            M[center, east] = 1
            M[center, center] = 1
        elseif i == N
            M[center, west] = 1
            M[center, center] = 1
        else
            M[center, east] = 1
            M[center, west] = 1
            M[center, center] = 2
        end
    end
    M .*= dx
    return M, b
end

function continuity_x(u, mesh::CartesianMesh)
    N = mesh.N
    dx = 1 / N
    r = zeros(N^2)
    for i in 1:N, j in 1:N
        if i == 1
            uw = zero(u[1])
            ue = u[cartesian_to_linear_indices(i + 1, j, mesh)]
        elseif i == N
            ue = zero(u[N])
            uw = u[cartesian_to_linear_indices(i - 1, j, mesh)]
        else
            uw = u[cartesian_to_linear_indices(i - 1, j, mesh)]
            ue = u[cartesian_to_linear_indices(i + 1, j, mesh)]
        end
        r[cartesian_to_linear_indices(i, j, mesh)] = (ue - uw)
    end
    r ./= dx
    return r
end

function continuity_y(v, mesh::CartesianMesh)
    N = mesh.N
    dy = 1 / N
    r = zeros(N^2)
    for i in 1:N, j in 1:N
        if j == 1
            vs = zero(v[1])
            vn = v[cartesian_to_linear_indices(i, j + 1, mesh)]
        elseif j == N
            vn = zero(v[N])
            vs = v[cartesian_to_linear_indices(i, j - 1, mesh)]
        else
            vs = v[cartesian_to_linear_indices(i, j - 1, mesh)]
            vn = v[cartesian_to_linear_indices(i, j + 1, mesh)]
        end
        r[cartesian_to_linear_indices(i, j, mesh)] = (vn - vs)
    end
    r ./= dy
    return r
end

function get_pressure_dx_operator(mesh::CartesianMesh, type::Type{T}=Float64) where T<:AbstractFloat
    N = mesh.N
    dx = 1 / N
    M = spzeros(T, N^2, N^2)
    b = spzeros(T, N^2)
    for i in 1:N, j in 1:N
        center = cartesian_to_linear_indices(i, j, mesh)
        east = i != N ? cartesian_to_linear_indices(i + 1, j, mesh) : nothing
        west = i != 1 ? cartesian_to_linear_indices(i - 1, j, mesh) : nothing
        if i == 1
            M[center, east] = 1
            M[center, center] -= 1
        elseif i == N
            M[center, west] = -1
            M[center, center] += 1
        else
            M[center, east] = 1
            M[center, west] = -1
        end
    end
    M ./= (2 * dx)
    return M, b
end

function pressure_dx(p, mesh::CartesianMesh)
    N = mesh.N
    dx = 1 / N
    r = zeros(N^2)
    for i in 1:N, j in 1:N
        center = cartesian_to_linear_indices(i, j, mesh)
        east = i != N ? cartesian_to_linear_indices(i + 1, j, mesh) : nothing
        west = i != 1 ? cartesian_to_linear_indices(i - 1, j, mesh) : nothing
        if i == 1
            dpw = zero(p[1])
            dpe = p[east] - p[center]
        elseif i == N
            dpw = p[center] - p[west]
            dpe = zero(p[N])
        else
            dpw = p[center] - p[west]
            dpe = p[east] - p[center]
        end
        r[center] = (dpe + dpw)
    end
    r ./= 2 * dx
    return r
end

function get_pressure_dy_operator(mesh::CartesianMesh, type::Type{T}=Float64) where T<:AbstractFloat
    N = mesh.N
    dy = 1 / N
    M = spzeros(T, N^2, N^2)
    b = spzeros(T, N^2)
    for j in 1:N, i in 1:N
        center = cartesian_to_linear_indices(i, j, mesh)
        north = j != N ? cartesian_to_linear_indices(i, j + 1, mesh) : nothing
        south = j != 1 ? cartesian_to_linear_indices(i, j - 1, mesh) : nothing
        if j == 1
            M[center, north] = 1
            M[center, center] -= 1
        elseif j == N
            M[center, south] = -1
            M[center, center] -= 1
        else
            M[center, north] = 1
            M[center, south] = -1
        end
    end
    M ./= (2 * dy)
    return M, b
end

function pressure_dy(p, mesh::CartesianMesh)
    N = mesh.N
    dy = 1 / N
    r = zeros(N^2)
    zero_val = zero(eltype(p))
    for i in 1:N, j in 1:N
        center = cartesian_to_linear_indices(i, j, mesh)
        north = cartesian_to_linear_indices(i, j + 1, mesh)
        south = cartesian_to_linear_indices(i, j - 1, mesh)
        if j == 1
            dps = zero_val
            dpn = p[north] - p[center]
        elseif j == N
            dps = p[center] - p[south]
            dpn = -2 * p[center]
        else
            dps = p[center] - p[south]
            dpn = p[north] - p[center]
        end
        r[center] = (dpn + dps)
    end
    r ./= 2 * dy
    return r
end

function get_laplacian_operator(mesh::CartesianMesh, topvalue::Number)
    """
    assemble laplacian operator Matrix and b vector such that Δg = Mg + b
    """
    N = mesh.N
    dx = 1 / N
    dy = 1 / N
    M = spzeros(Float64, N^2, N^2)
    b = spzeros(Float64, N^2)

    for j in 1:N, i in 1:N
        center = cartesian_to_linear_indices(i, j, mesh)
        j != N ? north = cartesian_to_linear_indices(i, j + 1, mesh) : north = nothing
        j != 1 ? south = cartesian_to_linear_indices(i, j - 1, mesh) : south = nothing
        i != N ? east = cartesian_to_linear_indices(i + 1, j, mesh) : east = nothing
        i != 1 ? west = cartesian_to_linear_indices(i - 1, j, mesh) : west = nothing
        if i == 1
            M[center, east] = 1
            M[center, center] -= 2
        elseif i == N
            M[center, west] = 1
            M[center, center] -= 2
        else
            M[center, west] = 1
            M[center, center] -= 2
            M[center, east] = 1
        end
        if j == 1
            M[center, north] = 1
            M[center, center] -= 2
        elseif j == N
            M[center, south] = 1
            M[center, center] -= 2
            b[center] = topvalue
        else
            M[center, north] = 1
            M[center, center] -= 2
            M[center, south] = 1
        end
    end
    M ./= dx * dy
    b ./= dx * dy
    return M, b
end

"""
    compute_face_values(ϕ, mesh, u, v, limiter, bc=DomainBoundaryConditions())

Compute TVD face values using multiple dispatch for clean, extensible code.
Supports different boundary conditions on different sides.
"""
function compute_face_values(ϕ, mesh::CartesianMesh, u, v, limiter::AbstractLimiter,
    bc::DomainBoundaryConditions=DomainBoundaryConditions())
    N = mesh.N
    total_faces = get_total_faces(mesh)
    ϕ_f = similar(ϕ, total_faces)

    @inbounds for face_idx in 1:total_faces
        ϕ_f[face_idx] = compute_single_face_value(ϕ, mesh, u, v, limiter, bc, face_idx)
    end

    return ϕ_f
end

function compute_single_face_value(ϕ, mesh::CartesianMesh, u, v, limiter::AbstractLimiter,
    bc::DomainBoundaryConditions, face_idx::Int)
    i1, j1, i2, j2, orientation = face_neighbors(mesh, face_idx)

    # Get cell indices and values
    cell1 = cartesian_to_linear_indices(i1, j1, mesh)
    cell2 = cartesian_to_linear_indices(i2, j2, mesh)

    # Determine face velocity and flow direction using dispatch
    face_orient = orientation == :vertical ? Vertical() : Horizontal()
    u_face = 0.5 * (face_velocity(u[cell1], v[cell1], face_orient) +
                    face_velocity(u[cell2], v[cell2], face_orient))

    dir = flow_direction(u_face)

    # Get upwind/downwind cells using dispatch
    upwind_idx = upwind_cell(cell1, cell2, dir)
    downwind_idx = downwind_cell(cell1, cell2, dir)
    upwind_i, upwind_j = upwind_idx == cell1 ? (i1, j1) : (i2, j2)
    downwind_i, downwind_j = downwind_idx == cell1 ? (i1, j1) : (i2, j2)

    # Compute TVD face value
    ϕ_upwind = ϕ[upwind_idx]
    ϕ_downwind = ϕ[downwind_idx]

    r = compute_r_ratio(ϕ, mesh, upwind_i, upwind_j, downwind_i, downwind_j,
        face_orient, bc)

    return ϕ_upwind + 0.5 * limiter(r) * (ϕ_downwind - ϕ_upwind)
end

"""
    compute_r_ratio(ϕ, mesh, upwind_i, upwind_j, downwind_i, downwind_j, orientation, bc)

Compute r-ratio using multiple dispatch for boundary conditions.
Supports different boundary conditions on different sides.
"""
function compute_r_ratio(ϕ, mesh::CartesianMesh, upwind_i, upwind_j, downwind_i, downwind_j,
    orientation::FaceOrientation, bc::DomainBoundaryConditions)
    N = mesh.N

    # Find upwind-upwind cell using dispatch
    upwind_upwind_i, upwind_upwind_j = get_upwind_upwind_cell(
        upwind_i, upwind_j, downwind_i, downwind_j, orientation, N)

    # Get values
    ϕ_upwind = ϕ[cartesian_to_linear_indices(upwind_i, upwind_j, mesh)]
    ϕ_downwind = ϕ[cartesian_to_linear_indices(downwind_i, downwind_j, mesh)]
    ϕ_upwind_upwind = get_upwind_upwind_value(ϕ, mesh, upwind_upwind_i, upwind_upwind_j,
        upwind_i, upwind_j, bc)

    # Compute r ratio with numerical stability
    denominator = ϕ_downwind - ϕ_upwind
    return abs(denominator) < 1e-12 ? 1.0 : (ϕ_upwind - ϕ_upwind_upwind) / denominator
end

# Multiple dispatch for different face orientations
function get_upwind_upwind_cell(upwind_i, upwind_j, downwind_i, downwind_j,
    ::Vertical, N)
    if downwind_i > upwind_i  # Flow left to right
        return max(1, upwind_i - 1), upwind_j
    else  # Flow right to left
        return min(N, upwind_i + 1), upwind_j
    end
end

function get_upwind_upwind_cell(upwind_i, upwind_j, downwind_i, downwind_j,
    ::Horizontal, N)
    if downwind_j > upwind_j  # Flow bottom to top
        return upwind_i, max(1, upwind_j - 1)
    else  # Flow top to bottom
        return upwind_i, min(N, upwind_j + 1)
    end
end

# Multiple dispatch for boundary conditions with domain-specific BCs
function get_upwind_upwind_value(ϕ, mesh, upwind_upwind_i, upwind_upwind_j,
    upwind_i, upwind_j, bc::DomainBoundaryConditions)
    N = mesh.N
    if upwind_upwind_i == upwind_i && upwind_upwind_j == upwind_j
        # At boundary - determine which side and get appropriate BC
        boundary_side = get_boundary_side(upwind_i, upwind_j, N)
        if boundary_side !== nothing
            specific_bc = get_bc(bc, boundary_side)
            return get_boundary_value(ϕ, mesh, upwind_i, upwind_j, specific_bc)
        else
            # Interior cell (shouldn't happen in this context)
            return ϕ[cartesian_to_linear_indices(upwind_i, upwind_j, mesh)]
        end
    else
        return ϕ[cartesian_to_linear_indices(upwind_upwind_i, upwind_upwind_j, mesh)]
    end
end

# Helper function to get boundary value based on BC type
function get_boundary_value(ϕ, mesh, i, j, bc::DirichletBC)
    return bc.value
end

function get_boundary_value(ϕ, mesh, i, j, bc::NeumannBC)
    # Zero gradient extrapolation
    return ϕ[cartesian_to_linear_indices(i, j, mesh)]
end

"""
get_convection_operators(mesh, ϕ, u, v, limiter, bc)

returns the convection operators for u and v using TVD scheme
Supports different boundary conditions on different sides.
"""
function set_convection_operators!(M_u, M_v, mesh::CartesianMesh, ϕ, u, v, limiter::AbstractLimiter,
    bc::DomainBoundaryConditions=DomainBoundaryConditions())
    N = mesh.N
    dx = dy = 1 / N

    ϕ_f = compute_face_values(ϕ, mesh, u, v, limiter, bc)

    for i in 1:N, j in 1:N
        center = cartesian_to_linear_indices(i, j, mesh)

        # Process each direction using helper functions
        process_convection_direction!(M_u, M_v, ϕ_f, mesh, i, j, center, dx, dy)
    end
    nothing
end

function process_convection_direction!(M_u, M_v, ϕ_f, mesh, i, j, center, dx, dy)
    N = mesh.N

    # East/West faces for u-momentum
    if i < N
        east_face_idx = get_face_index(mesh, i, j, :east)
        east_neighbor = cartesian_to_linear_indices(i + 1, j, mesh)
        ϕ_e = ϕ_f[east_face_idx]
        M_u[center, east_neighbor] += ϕ_e * dx / 2
        M_u[center, center] += ϕ_e * dx / 2
    end

    if i > 1
        west_face_idx = get_face_index(mesh, i, j, :west)
        west_neighbor = cartesian_to_linear_indices(i - 1, j, mesh)
        ϕ_w = ϕ_f[west_face_idx]
        M_u[center, west_neighbor] -= ϕ_w * dx / 2
        M_u[center, center] -= ϕ_w * dx / 2
    end

    # North/South faces for v-momentum
    if j < N
        north_face_idx = get_face_index(mesh, i, j, :north)
        north_neighbor = cartesian_to_linear_indices(i, j + 1, mesh)
        ϕ_n = ϕ_f[north_face_idx]
        M_v[center, north_neighbor] += ϕ_n * dy / 2
        M_v[center, center] += ϕ_n * dy / 2
    end

    if j > 1
        south_face_idx = get_face_index(mesh, i, j, :south)
        south_neighbor = cartesian_to_linear_indices(i, j - 1, mesh)
        ϕ_s = ϕ_f[south_face_idx]
        M_v[center, south_neighbor] -= ϕ_s * dy / 2
        M_v[center, center] -= ϕ_s * dy / 2
    end
end