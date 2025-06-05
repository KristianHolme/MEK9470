function split_uvp(uvp, mesh::CartesianMesh)
    N = mesh.N
    @assert length(uvp) == 3 * N^2
    u = uvp[1:N^2]
    v = uvp[N^2+1:2*N^2]
    p = uvp[2*N^2+1:3*N^2]
    return u, v, p
end

"""
continuity in x direction(for u).

applying the divergence theorem, the continuity equation is requiring the sum of the fluxes out of each cell to be zero.

here we make the x-direction contribution to the continuity equation, and the y-direction contribution is the transpose of this.

at the left side (i==1) the velocity at the the left side is zero and we approximate the flux at the boundary using interpolation, 
resulting in the face flux being the average of u_E and u_P. so the x-part of the equation becomes (u_E + u_P)/2 *dy

for i == N, on the right side, the flux u_e = 0, so we only get the left term -(u_W + u_P)/2 *dy

internally we get the sum of the previous terms, equaling (u_E + u_P)/2 *dy - (u_W + u_P)/2 *dy = (u_E - u_W)/2 *dy
"""
function get_continuity_operator_x(mesh::CartesianMesh, type::Type{T}=Float64) where T<:AbstractFloat
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
            M[center, west] = -1
            M[center, center] = -1
        else
            M[center, east] = 1
            M[center, west] = -1
            # M[center, center] = 2
        end
    end
    M .*= dx / 2
    return M, b
end


function get_continuity_operator_y(mesh::CartesianMesh, type::Type{T}=Float64) where T<:AbstractFloat
    N = mesh.N
    dy = 1 / N
    M = spzeros(T, N^2, N^2)
    b = spzeros(T, N^2)
    for i in 1:N, j in 1:N
        center = cartesian_to_linear_indices(i, j, mesh)
        north = j != N ? cartesian_to_linear_indices(i, j + 1, mesh) : nothing
        south = j != 1 ? cartesian_to_linear_indices(i, j - 1, mesh) : nothing

        if j == 1
            M[center, north] = 1
            M[center, center] = 1
        elseif j == N
            M[center, south] = -1
            M[center, center] = -1
        else
            M[center, north] = 1
            M[center, south] = -1
            # M[center, center] = 2
        end
    end
    M .*= dy / 2
    return M, b
end

"""
get the pressure gradient operator in the x-direction.
we have neumann bc on all sides, ∂p/∂x = 0

we calculate this as te area of the cell (Δx * Δy) times the difference in pressure at the cell faces.

the pressure gradient of the cell is computed as the average of the interpolated pressure at the east and west faces:
∫∂p/∂x dxdy ≈ ΔxΔy(∂p/∂x|_e + ∂p/∂x |_w)/2 ≈ ΔxΔy( (p_E - p_P)/Δx + (p_P - p_W)/Δx)/2 ≈ ΔxΔy(p_E - p_W)/2Δx = Δy(p_E - p_W)/2

on the boundaries, one of the terms is zero, so for the east side we get ΔxΔy( 0 + (p_P - p_W)/Δx)/2 = Δy(p_P - p_W)/2,
and for the west side we get ΔxΔy( (p_E - p_P)/Δx + 0)/2 = Δy(p_E - p_P)/2
"""
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
            M[center, center] = -1
        elseif i == N
            M[center, center] = 1
            M[center, west] = -1
        else
            M[center, east] = 1
            M[center, west] = -1
        end
    end
    M .*= dx / 2
    return M, b
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
            M[center, center] = -1
        elseif j == N
            M[center, center] = 1
            M[center, south] = -1
        else
            M[center, north] = 1
            M[center, south] = -1
        end
    end
    M .*= dy / 2
    return M, b
end

"""
assemble laplacian operator Matrix and b vector such that Δg = Mg + b

in the x-direction we want the second derivative. We express this as the finite difference of the finite difference 
derivatives at the cell faces: Δu|P ≈ (∂u/∂x|e - ∂u/∂x|w)/dx ≈ [(u_E - u_P)/dx - (u_P - u_W)/dx]/dx = (u_E - 2u_P + u_W)/dx^2

on the left boundary, we have u_w = 0, so ∂u/∂x|w ≈ (u_P - u_w)/(dx/2) = 2u_P/dx, and we get the term 
    (∂u/∂x|e - ∂u/∂x|w)/dx ≈ [(u_E - u_P)/dx - 2u_P/dx]/dx = (u_E - 3u_P)/dx^2

on the right boundary, we have u_e = 0, so ∂u/∂x|e ≈ (u_e - u_P)/(dx/2) = -2u_P/dx, and we get the term 
    (∂u/∂x|e - ∂u/∂x|w)/dx ≈ [-2u_P/dx - (u_P - u_W)/dx]/dx = (-3u_P + u_W)/dx^2

also when integrating over a cell, we multiply by the area of the cell (Δx * Δy) = dx^2
"""
function get_laplacian_operator(mesh::CartesianMesh, topvalue::Number)
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

        # X-direction terms
        if i == 1
            M[center, east] = 1
            M[center, center] -= 3
        elseif i == N
            M[center, west] = 1
            M[center, center] -= 3
        else
            M[center, west] = 1
            M[center, center] -= 2
            M[center, east] = 1
        end

        # Y-direction terms
        if j == 1
            M[center, north] = 1
            M[center, center] -= 3
        elseif j == N
            M[center, south] = 1
            M[center, center] -= 3
            ccx, ccy = cell_center(i, j, mesh)
            regular_topvalue = 16 * (1 - ccx)^2 * ccx^2 * topvalue
            b[center] = regular_topvalue * 2
        else
            M[center, north] = 1
            M[center, center] -= 2
            M[center, south] = 1
        end
    end
    # M ./= dx * dy
    # b ./= dx * dy
    return M, b
end

"""
    compute_face_values(ϕ, mesh, u, v, limiter, bc=DomainBoundaryConditions())

Compute TVD face values.
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

    any(isnan.(ϕ_f)) && @warn "ϕ_f contains nans"

    return ϕ_f
end

function compute_single_face_value(ϕ, mesh::CartesianMesh, u, v, limiter::AbstractLimiter,
    bc::DomainBoundaryConditions, face_idx::Int)
    i1, j1, i2, j2, orientation = face_neighbors(mesh, face_idx)

    # Get cell indices and values
    cell1 = cartesian_to_linear_indices(i1, j1, mesh)
    cell2 = cartesian_to_linear_indices(i2, j2, mesh)

    # Determine face velocity
    face_orient = orientation == :vertical ? Vertical() : Horizontal()
    u_face = 0.5 * (face_velocity(u[cell1], v[cell1], face_orient) +
                    face_velocity(u[cell2], v[cell2], face_orient))

    # Deterministic flow direction selection
    if u_face >= 0
        # Positive flow: cell1 -> cell2 (upwind=cell1, downwind=cell2)
        upwind_i, upwind_j = i1, j1
        downwind_i, downwind_j = i2, j2
        ϕ_upwind = ϕ[cell1]
        ϕ_downwind = ϕ[cell2]
    else
        # Negative flow: cell2 -> cell1 (upwind=cell2, downwind=cell1)
        upwind_i, upwind_j = i2, j2
        downwind_i, downwind_j = i1, j1
        ϕ_upwind = ϕ[cell2]
        ϕ_downwind = ϕ[cell1]
    end

    # Compute r-ratio for the determined flow direction
    r = compute_r_ratio(ϕ, mesh, upwind_i, upwind_j, downwind_i, downwind_j, face_orient, bc)

    # Compute TVD face value
    ϕ_face = ϕ_upwind + 0.5 * limiter(r) * (ϕ_downwind - ϕ_upwind)

    return ϕ_face
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
    denominator = max(denominator, 1e-8)
    r = (ϕ_upwind - ϕ_upwind_upwind) / denominator
    return r
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
            return get_upwind_upwind_bc(ϕ, mesh, upwind_i, upwind_j, specific_bc)
        else
            # Interior cell (shouldn't happen in this context)
            return ϕ[cartesian_to_linear_indices(upwind_i, upwind_j, mesh)]
        end
    else
        return ϕ[cartesian_to_linear_indices(upwind_upwind_i, upwind_upwind_j, mesh)]
    end
end

# Helper function to get boundary value based on BC type
function get_upwind_upwind_bc(ϕ, mesh, i, j, bc::DirichletBC)
    #using ghost cell extrapolation
    boundary_cell_value = ϕ[cartesian_to_linear_indices(i, j, mesh)]
    if bc.value > 0.5
        cxx, _ = cell_center(i, j, mesh)
        regular_topvalue = 16 * (1 - cxx)^2 * cxx^2
        bc_value = regular_topvalue
    else
        bc_value = bc.value
    end
    upwind_upwind_value = 2 * bc_value - boundary_cell_value
    return upwind_upwind_value
end

function get_upwind_upwind_bc(ϕ, mesh, i, j, bc::NeumannBC)
    # Zero gradient extrapolation
    return ϕ[cartesian_to_linear_indices(i, j, mesh)]
end

"""
get_convection_operators(mesh, ϕ, u, v, limiter, bc)

returns the convection operators for u and v using TVD scheme
Supports different boundary conditions on different sides.
"""
function set_convection_operators!(M_u, M_v, mesh::CartesianMesh, ϕ_f, u, v, limiter::AbstractLimiter,
    bc::DomainBoundaryConditions=DomainBoundaryConditions())
    N = mesh.N
    dx = dy = 1 / N

    Iu = Int[]
    Ju = Int[]
    Vu = Float64[]
    Iv = Int[]
    Jv = Int[]
    Vv = Float64[]


    for i in 1:N, j in 1:N
        # Process each direction using helper functions
        process_convection_direction!(Iu, Ju, Vu, Iv, Jv, Vv, ϕ_f, mesh, i, j)
    end
    M_u = sparse(Iu, Ju, Vu .* dx * 0.5, N^2, N^2)
    M_v = sparse(Iv, Jv, Vv .* dy * 0.5, N^2, N^2)
    # M_u .*= dx * 0.5
    # M_v .*= dy * 0.5
    nothing
end

function process_convection_direction!(Iu, Ju, Vu, Iv, Jv, Vv, ϕ_f, mesh, i, j)
    N = mesh.N
    center = cartesian_to_linear_indices(i, j, mesh)

    # East/West faces for u-momentum
    center_value_u = 0.0
    center_value_v = 0.0
    if i < N
        east_face_idx = get_face_index(mesh, i, j, :east)
        east_neighbor = cartesian_to_linear_indices(i + 1, j, mesh)
        ϕ_e = ϕ_f[east_face_idx]
        # M_u[center, east_neighbor] += ϕ_e
        # M_u[center, center] += ϕ_e
        push!(Iu, center)
        push!(Ju, east_neighbor)
        push!(Vu, ϕ_e)
        center_value_u += ϕ_e
    end

    if i > 1
        west_face_idx = get_face_index(mesh, i, j, :west)
        west_neighbor = cartesian_to_linear_indices(i - 1, j, mesh)
        ϕ_w = ϕ_f[west_face_idx]
        # M_u[center, west_neighbor] -= ϕ_w
        # M_u[center, center] -= ϕ_w
        push!(Iu, center)
        push!(Ju, west_neighbor)
        push!(Vu, -ϕ_w)
        center_value_u -= ϕ_w
    end

    # North/South faces for v-momentum
    if j < N
        north_face_idx = get_face_index(mesh, i, j, :north)
        north_neighbor = cartesian_to_linear_indices(i, j + 1, mesh)
        ϕ_n = ϕ_f[north_face_idx]
        # M_v[center, north_neighbor] += ϕ_n
        # M_v[center, center] += ϕ_n
        push!(Iv, center)
        push!(Jv, north_neighbor)
        push!(Vv, ϕ_n)
        center_value_v += ϕ_n
    end

    if j > 1
        south_face_idx = get_face_index(mesh, i, j, :south)
        south_neighbor = cartesian_to_linear_indices(i, j - 1, mesh)
        ϕ_s = ϕ_f[south_face_idx]
        # M_v[center, south_neighbor] -= ϕ_s
        # M_v[center, center] -= ϕ_s
        push!(Iv, center)
        push!(Jv, south_neighbor)
        push!(Vv, -ϕ_s)
        center_value_v -= ϕ_s
    end
    push!(Iu, center)
    push!(Ju, center)
    push!(Vu, center_value_u)
    push!(Iv, center)
    push!(Jv, center)
    push!(Vv, center_value_v)
    nothing
end

function get_convection_correction_values(mwi_face_corrections, ϕu_f, ϕv_f, mesh::CartesianMesh)
    N = mesh.N
    b_u_corrections = zeros(eltype(mwi_face_corrections), N^2)
    b_v_corrections = zeros(eltype(mwi_face_corrections), N^2)
    set_convection_corrections!(b_u_corrections, ϕu_f, mwi_face_corrections, mesh)
    set_convection_corrections!(b_v_corrections, ϕv_f, mwi_face_corrections, mesh)
    return b_u_corrections, b_v_corrections
end

function set_convection_corrections!(b_corrections, ϕ_f, mwi_face_corrections, mesh::CartesianMesh)
    N = mesh.N
    dx = dy = 1 / N
    for i in 1:N, j in 1:N
        center = cartesian_to_linear_indices(i, j, mesh)

        # East/West faces for u-momentum
        if i < N
            east_face_idx = get_face_index(mesh, i, j, :east)
            east_neighbor = cartesian_to_linear_indices(i + 1, j, mesh)
            ϕ_e = ϕ_f[east_face_idx]
            b_corrections[center] += ϕ_e * dx * mwi_face_corrections[east_face_idx]
        end

        if i > 1
            west_face_idx = get_face_index(mesh, i, j, :west)
            west_neighbor = cartesian_to_linear_indices(i - 1, j, mesh)
            ϕ_w = ϕ_f[west_face_idx]
            b_corrections[center] -= ϕ_w * dx * mwi_face_corrections[west_face_idx]
        end

        # North/South faces for v-momentum
        if j < N
            north_face_idx = get_face_index(mesh, i, j, :north)
            north_neighbor = cartesian_to_linear_indices(i, j + 1, mesh)
            ϕ_n = ϕ_f[north_face_idx]
            b_corrections[center] += ϕ_n * dy * mwi_face_corrections[north_face_idx]
        end

        if j > 1
            south_face_idx = get_face_index(mesh, i, j, :south)
            south_neighbor = cartesian_to_linear_indices(i, j - 1, mesh)
            ϕ_s = ϕ_f[south_face_idx]
            b_corrections[center] -= ϕ_s * dy * mwi_face_corrections[south_face_idx]
        end
    end
end
"""
get the Rhie & Chow correction values for each face.

U_f = ̄U_f - ̄d_f(∂p/∂x|_f i \bar{∂p/∂x}|_f)

this function returns d_f(∂p/∂x|_f i \bar{∂p/∂x}|_f) for each face (derivative in y direction for horizontal faces, x direction for vertical faces)
"""
function get_face_correction_values(mesh::CartesianMesh, uvp::AbstractVector, momentum_u, momentum_v, p_dx, p_dy)
    N = mesh.N
    dx = dy = 1 / N
    u, v, p = split_uvp(uvp, mesh)
    total_faces = get_total_faces(mesh)
    corrections = zeros(eltype(uvp), total_faces)

    # Pressure gradients (multiplied by cell volume)
    pressure_gradient_x = p_dx * p
    pressure_gradient_y = p_dy * p

    for f in 1:total_faces
        i1, j1, i2, j2, orientation = face_neighbors(mesh, f)
        neighbor1 = cartesian_to_linear_indices(i1, j1, mesh)
        neighbor2 = cartesian_to_linear_indices(i2, j2, mesh)
        # face gradient (P_F - P_P) / dx multiplied by volume (dx^2)
        dp_f = (p[neighbor2] - p[neighbor1]) * dx
        if orientation == :Vertical
            # flux is in the x direction
            a_P = momentum_u[neighbor1, neighbor1]
            a_F = momentum_u[neighbor2, neighbor2]
            dp_f_1 = pressure_gradient_x[neighbor1]
            dp_f_2 = pressure_gradient_x[neighbor2]
        else
            # flux is in the y direction
            a_P = momentum_v[neighbor1, neighbor1]
            a_F = momentum_v[neighbor2, neighbor2]
            dp_f_1 = pressure_gradient_y[neighbor1]
            dp_f_2 = pressure_gradient_y[neighbor2]
        end
        dP = 1 / a_P
        dF = 1 / a_F
        df = 0.5 * (dP + dF)
        dp_f_interp = (dp_f_1 + dp_f_2) * 0.5
        # a_f = (a_P + a_F) / 2

        # df = 1 / a_f
        correction = df * (dp_f - dp_f_interp)
        corrections[f] = correction
    end
    return corrections
end



function get_continuity_correction_values(mesh::CartesianMesh, uvp, Muu, Mvv, p_dx, p_dy)
    face_corrections = get_face_correction_values(mesh, uvp, Muu, Mvv, p_dx, p_dy)
    return get_continuity_correction_values(face_corrections, mesh)
end

function get_continuity_correction_values(face_corrections, mesh::CartesianMesh)
    N = mesh.N
    dx = dy = 1 / N
    b = zeros(eltype(face_corrections), mesh.N^2)
    # any(isnan.(face_corrections)) && @warn "face_corrections contains nans"
    for i in 1:mesh.N, j in 1:mesh.N
        center = cartesian_to_linear_indices(i, j, mesh)
        north_face_idx = j != N ? get_face_index(mesh, i, j, :north) : nothing
        south_face_idx = j != 1 ? get_face_index(mesh, i, j, :south) : nothing
        east_face_idx = i != N ? get_face_index(mesh, i, j, :east) : nothing
        west_face_idx = i != 1 ? get_face_index(mesh, i, j, :west) : nothing
        if i == 1
            b[center] += face_corrections[east_face_idx]
        elseif i == mesh.N
            b[center] -= face_corrections[west_face_idx]
        else
            b[center] += face_corrections[east_face_idx] - face_corrections[west_face_idx]
        end
        if j == 1
            b[center] += face_corrections[north_face_idx]
        elseif j == mesh.N
            b[center] -= face_corrections[south_face_idx]
        else
            b[center] += face_corrections[north_face_idx] - face_corrections[south_face_idx]
        end
    end
    # Corrected scaling: multiply by face area (dx for the face normal to flux direction)
    b .*= dx
    return b
end