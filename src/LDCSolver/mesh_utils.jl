function linear_to_cartesian_indices(ind::Int, mesh::CartesianMesh)
    N = mesh.N
    i = (ind - 1) % N + 1
    j = (ind - 1) ÷ N + 1
    @assert i <= N && j <= N "Index out of bounds"
    return i, j
end

function cartesian_to_linear_indices(i::Int, j::Int, mesh::CartesianMesh)
    N = mesh.N
    @assert i <= N && j <= N "Index out of bounds"
    return (j - 1) * N + i
end

function cell_center(i::Int, j::Int, mesh::CartesianMesh)
    N = mesh.N
    dx = 1/N
    dy = 1/N
    x = (i - 1)*dx + dx/2
    y = (j - 1)*dy + dy/2
    return (x, y)
end

cell_center(ind::Int, mesh::CartesianMesh) = cell_center(linear_to_cartesian_indices(ind, mesh)..., mesh)

