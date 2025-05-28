abstract type AbstractMesh end

struct CartesianMesh <: AbstractMesh
    N::Int
end

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
    dx = 1 / N
    dy = 1 / N
    x = (i - 1) * dx + dx / 2
    y = (j - 1) * dy + dy / 2
    return (x, y)
end

cell_center(ind::Int, mesh::CartesianMesh) = cell_center(linear_to_cartesian_indices(ind, mesh)..., mesh)

"""
    get_face_index(mesh::CartesianMesh, i::Int, j::Int, direction::Symbol)

Get the linear index of a face for cell (i,j) in the given direction.

Face indexing scheme:
- Vertical faces (east/west): indexed 1 to N*(N-1)
- Horizontal faces (north/south): indexed N*(N-1)+1 to 2*N*(N-1)

For vertical faces:
- East face of cell (i,j): face between cells (i,j) and (i+1,j)
- West face of cell (i,j): face between cells (i-1,j) and (i,j)

For horizontal faces:
- North face of cell (i,j): face between cells (i,j) and (i,j+1)  
- South face of cell (i,j): face between cells (i,j-1) and (i,j)

Arguments:
- mesh: CartesianMesh
- i, j: cell coordinates (1-indexed)
- direction: :east, :west, :north, or :south

Returns:
- Linear face index
"""
function get_face_index(mesh::CartesianMesh, i::Int, j::Int, direction::Symbol)
    N = mesh.N
    @assert 1 <= i <= N && 1 <= j <= N "Cell indices out of bounds"

    if direction == :east
        @assert i < N "No east face for rightmost cells"
        # East face between cells (i,j) and (i+1,j)
        # Vertical face at x-position i+0.5, y-position j
        return (j - 1) * (N - 1) + i

    elseif direction == :west
        @assert i > 1 "No west face for leftmost cells"
        # West face between cells (i-1,j) and (i,j)
        # This is the same as the east face of cell (i-1,j)
        return (j - 1) * (N - 1) + (i - 1)

    elseif direction == :north
        @assert j < N "No north face for topmost cells"
        # North face between cells (i,j) and (i,j+1)
        # Horizontal face at x-position i, y-position j+0.5
        return N * (N - 1) + (j - 1) * N + i

    elseif direction == :south
        @assert j > 1 "No south face for bottommost cells"
        # South face between cells (i,j-1) and (i,j)
        # This is the same as the north face of cell (i,j-1)
        return N * (N - 1) + (j - 2) * N + i

    else
        error("Invalid direction: $direction. Must be :east, :west, :north, or :south")
    end
end

"""
    get_total_faces(mesh::CartesianMesh)

Get the total number of internal faces in the mesh.

For an NxN grid:
- Vertical faces: N*(N-1)
- Horizontal faces: (N-1)*N  
- Total: 2*N*(N-1)
"""
function get_total_faces(mesh::CartesianMesh)
    N = mesh.N
    return 2 * N * (N - 1)
end

"""
    face_neighbors(mesh::CartesianMesh, face_idx::Int)

Get the two cells that share a face.

Returns:
- (i1, j1, i2, j2): coordinates of the two neighboring cells
- direction: :vertical or :horizontal
"""
function face_neighbors(mesh::CartesianMesh, face_idx::Int)
    N = mesh.N
    total_vertical_faces = N * (N - 1)

    if face_idx <= total_vertical_faces
        # Vertical face
        j = div(face_idx - 1, N - 1) + 1
        i = mod(face_idx - 1, N - 1) + 1
        return (i, j, i + 1, j, :vertical)
    else
        # Horizontal face
        horizontal_idx = face_idx - total_vertical_faces
        j = div(horizontal_idx - 1, N) + 1
        i = mod(horizontal_idx - 1, N) + 1
        return (i, j, i, j + 1, :horizontal)
    end
end


# Multiple dispatch for different boundary conditions
abstract type BoundaryCondition end
struct DirichletBC{T} <: BoundaryCondition
    value::T
end
struct NeumannBC{T} <: BoundaryCondition
    gradient::T
end

# Boundary side types for dispatch
abstract type BoundarySide end
struct North <: BoundarySide end
struct South <: BoundarySide end
struct East <: BoundarySide end
struct West <: BoundarySide end

# Comprehensive boundary condition specification for all sides
struct DomainBoundaryConditions{T}
    north::BoundaryCondition
    south::BoundaryCondition
    east::BoundaryCondition
    west::BoundaryCondition

    function DomainBoundaryConditions(north::BoundaryCondition, south::BoundaryCondition,
        east::BoundaryCondition, west::BoundaryCondition)
        return new{promote_type(typeof(north), typeof(south), typeof(east), typeof(west))}(
            north, south, east, west)
    end
end

# Convenience constructors for common cases
function DomainBoundaryConditions(; north=DirichletBC(0.0), south=DirichletBC(0.0),
    east=DirichletBC(0.0), west=DirichletBC(0.0))
    return DomainBoundaryConditions(north, south, east, west)
end

# Lid-driven cavity boundary conditions for u-velocity
function lid_driven_cavity_u_bc(lid_velocity=1.0)
    return DomainBoundaryConditions(
        north=DirichletBC(lid_velocity),  # Moving lid
        south=DirichletBC(0.0),           # No-slip bottom
        east=DirichletBC(0.0),            # No-slip right
        west=DirichletBC(0.0)             # No-slip left
    )
end

# Lid-driven cavity boundary conditions for v-velocity
function lid_driven_cavity_v_bc()
    return DomainBoundaryConditions(
        north=DirichletBC(0.0),   # No penetration at lid
        south=DirichletBC(0.0),   # No penetration at bottom
        east=DirichletBC(0.0),    # No penetration at right
        west=DirichletBC(0.0)     # No penetration at left
    )
end

# Get boundary condition for specific side using dispatch
get_bc(bc::DomainBoundaryConditions, ::North) = bc.north
get_bc(bc::DomainBoundaryConditions, ::South) = bc.south
get_bc(bc::DomainBoundaryConditions, ::East) = bc.east
get_bc(bc::DomainBoundaryConditions, ::West) = bc.west

# Determine which boundary side a cell is on
function get_boundary_side(i, j, N)
    if j == N
        return North()
    elseif j == 1
        return South()
    elseif i == N
        return East()
    elseif i == 1
        return West()
    else
        return nothing  # Interior cell
    end
end

# Face orientation types for dispatch
abstract type FaceOrientation end
struct Vertical <: FaceOrientation end
struct Horizontal <: FaceOrientation end

# Velocity component selection via dispatch
face_velocity(u, v, ::Vertical) = u
face_velocity(u, v, ::Horizontal) = v

# Direction types for cleaner dispatch
abstract type FlowDirection end
struct Positive <: FlowDirection end
struct Negative <: FlowDirection end

flow_direction(velocity) = velocity ≥ 0 ? Positive() : Negative()

# Upwind cell selection via dispatch
upwind_cell(cell1, cell2, ::Positive) = cell1
upwind_cell(cell1, cell2, ::Negative) = cell2
downwind_cell(cell1, cell2, dir::FlowDirection) = upwind_cell(cell2, cell1, dir)