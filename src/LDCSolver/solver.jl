function split_uvp(uvp, mesh::CartesianMesh)
    N = mesh.N
    @assert length(uvp) == 3*N^2
    u = uvp[1:N^2]
    v = uvp[N^2+1:2*N^2]
    p = uvp[2*N^2+1:3*N^2]
    return u, v, p
end

function continuity_equation(u, v, mesh::CartesianMesh)
    N = mesh.N
    dx = 1/N
    r = zeros(N^2)
    r .= continuity_x(u, mesh) + continuity_y(v, mesh)
    r .*= dx
    return r
end
"""
continuity in x direction(for u). Transpose for v in y-direction
"""
function get_continuity_operator(mesh::CartesianMesh, type::Type{T}=Float64) where T <: AbstractFloat
    N = mesh.N
    dx = 1/N
    M = spzeros(T,N^2, N^2)
    b = spzeros(T,N^2)
    for i in 1:N, j in 1:N
        center = cartesian_to_linear_indices(i, j, mesh)
        east = i != N ? cartesian_to_linear_indices(i+1, j, mesh) : nothing
        west = i != 1 ? cartesian_to_linear_indices(i+1, j, mesh) : nothing

        if i == 1
            M[center, east] = 1
        elseif i == N
            M[center, west] = -1
        else
            M[center, east] = 1
            M[center, west] = 1
        end
    end
    M ./= dx
    return M, b
end

function continuity_x(u, mesh::CartesianMesh)
    N = mesh.N
    r = zeros(N^2)
    for i in 1:N, j in 1:N
        if i == 1
            uw = zero(u[1])
            ue = u[cartesian_to_linear_indices(i+1, j, mesh)]
        elseif i == N
            ue = zero(u[N])
            uw = u[cartesian_to_linear_indices(i-1, j, mesh)]
        else
            uw = u[cartesian_to_linear_indices(i-1, j, mesh)]
            ue = u[cartesian_to_linear_indices(i+1, j, mesh)]
        end
        r[cartesian_to_linear_indices(i, j, mesh)] = (ue - uw)
    end
    r ./= dx
    return r
end

function continuity_y(v, mesh::CartesianMesh)
    N = mesh.N
    r = zeros(N^2)
    for i in 1:N, j in 1:N
        if j == 1
            vs = zero(v[1])
            vs = v[cartesian_to_linear_indices(i, j-1, mesh)]
        elseif j == N
            vn = zero(v[N])
            vn = v[cartesian_to_linear_indices(i, j+1, mesh)]
        else
            vs = v[cartesian_to_linear_indices(i, j-1, mesh)]
            vn = v[cartesian_to_linear_indices(i, j+1, mesh)]
        end
        r[cartesian_to_linear_indices(i, j, mesh)] = (vn -vs)
    end
    r ./= dy
    return r
end

function get_pressure_dx_operator(mesh::CartesianMesh, type::Type{T}=Float64) where T <: AbstractFloat
    N = mesh.N
    dx = 1/N
    M = spzeros(T,N^2, N^2)
    b = spzeros(T,N^2)
    for i in 1:N, j in 1:N
        center = cartesian_to_linear_indices(i, j, mesh)
        east = i != N ? cartesian_to_linear_indices(i+1, j, mesh) : nothing
        west = i != 1 ? cartesian_to_linear_indices(i-1, j, mesh) : nothing
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
    M ./= (2*dx)
    return M, b
end

function pressure_dx(p, mesh::CartesianMesh)
    N = mesh.N
    dx = 1/N
    r = zeros(N^2)
    for i in 1:N, j in 1:N
        center = cartesian_to_linear_indices(i, j, mesh)
        east = i != N ? cartesian_to_linear_indices(i+1, j, mesh) : nothing
        west = i != 1 ? cartesian_to_linear_indices(i-1, j, mesh) : nothing
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
    r ./= 2*dx
    return r
end

function get_pressure_dy_operator(mesh::CartesianMesh, type::Type{T}=Float64) where T <: AbstractFloat
    N = mesh.N
    dy = 1/N
    M = spzeros(T,N^2, N^2)
    b = spzeros(T,N^2)
    for j in 1:N, i in 1:N
        center = cartesian_to_linear_indices(i, j, mesh)
        north = j != N ? cartesian_to_linear_indices(i, j+1, mesh) : nothing
        south = j != 1 ? cartesian_to_linear_indices(i, j-1, mesh) : nothing
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
    M ./= (2*dy)
    return M, b
end

function pressure_dy(p, mesh::CartesianMesh)
    N = mesh.N
    dy = 1/N
    r = zeros(N^2)
    zero_val = zero(eltype(p))
    for i in 1:N, j in 1:N
        center = cartesian_to_linear_indices(i, j, mesh)
        north = cartesian_to_linear_indices(i, j+1, mesh)
        south = cartesian_to_linear_indices(i, j-1, mesh)
        if j == 1
            dps = zero_val
            dpn = p[north] - p[center]
        elseif j == N
            dps = p[center] - p[south]
            dpn = -2*p[center]
        else
            dps = p[center] - p[south]
            dpn = p[north] - p[center]
        end
        r[center] = (dpn + dps)
    end
    r ./= 2*dy
    return r
end

function get_laplacian_operator(mesh::CartesianMesh, topvalue::Number)
    """
    assemble laplacian operator Matrix and b vector such that Δg = Mg + b
    """
    N = mesh.N
    dx = 1/N
    dy = 1/N
    # r = zeros(N^2)
    M = spzeros(Float64,N^2, N^2)
    b = spzeros(Float64,N^2)

    for j in 1:N, i in 1:N
        center = cartesian_to_linear_indices(i, j, mesh)
        j != N ? north = cartesian_to_linear_indices(i, j+1, mesh) : north = nothing
        j != 1 ? south = cartesian_to_linear_indices(i, j-1, mesh) : south = nothing
        i != N ? east = cartesian_to_linear_indices(i+1, j, mesh) : east = nothing
        i != 1 ? west = cartesian_to_linear_indices(i-1, j, mesh) : west = nothing
        if i == 1
            # dxgw = g[center] # - zero
            # dxge = g[east] - g[center]
            M[center, east] = 1
            M[center, center] -= 2
        elseif i == N
            # dxgw = g[center] - g[west]
            # dxge = -g[center]
            M[center, west] = 1
            M[center, center] -= 2
        else
            # dxgw = g[center] - g[west]
            # dxge = g[east] - g[center]
            M[center, west] = 1
            M[center, center] -= 2
            M[center, east] = 1
        end
        if j == 1
            # dygn = g[north] - g[center]
            # dygs = g[center]
            M[center, north] = 1
            M[center, center] -= 2
        elseif j == N
            # dygn = topvalue - g[center]
            # dygs = g[center] - g[south]
            M[center, south] = 1
            M[center, center] -= 2
            b[center] = topvalue
        else
            # dygn = g[north] - g[center]
            # dygs = g[center] - g[south]
            M[center, north] = 1
            M[center, center] -= 2
            M[center, south] = 1
        end
        r[center] = (dxge - dxgw) + (dygn - dygs)
    end
    # r ./= dx*dy
    M ./= dx*dy
    b ./= dx*dy
    return M, b
end