using MEK9470

# Demonstrate the face indexing system for TVD schemes
println("=== Face Indexing System Demo ===")

# Create a small 3x3 mesh for demonstration
mesh = CartesianMesh(3)
N = 3

println("Mesh: $(N)x$(N) cells")
println("Total faces: $(get_total_faces(mesh))")
println("Vertical faces: $(N*(N-1)) = $(N*(N-1))")
println("Horizontal faces: $((N-1)*N) = $((N-1)*N)")
println()

# Show face indexing for each cell
println("=== Face Indices for Each Cell ===")
for j in N:-1:1  # Print from top to bottom
    for i in 1:N
        println("Cell ($i,$j):")

        # East face
        if i < N
            east_idx = get_face_index(mesh, i, j, :east)
            println("  East face: $east_idx")
        else
            println("  East face: boundary")
        end

        # West face  
        if i > 1
            west_idx = get_face_index(mesh, i, j, :west)
            println("  West face: $west_idx")
        else
            println("  West face: boundary")
        end

        # North face
        if j < N
            north_idx = get_face_index(mesh, i, j, :north)
            println("  North face: $north_idx")
        else
            println("  North face: boundary")
        end

        # South face
        if j > 1
            south_idx = get_face_index(mesh, i, j, :south)
            println("  South face: $south_idx")
        else
            println("  South face: boundary")
        end
        println()
    end
end

# Show face neighbors
println("=== Face Neighbors ===")
total_faces = get_total_faces(mesh)
for face_idx in 1:total_faces
    i1, j1, i2, j2, orientation = face_neighbors(mesh, face_idx)
    println("Face $face_idx ($orientation): connects cells ($i1,$j1) and ($i2,$j2)")
end
println()

# Demonstrate TVD face value computation
println("=== TVD Face Values Demo ===")

# Create a simple test case
ϕ = zeros(N^2)
u = zeros(N^2)
v = zeros(N^2)

# Set up a linear field in x-direction
for i in 1:N, j in 1:N
    idx = cartesian_to_linear_indices(i, j, mesh)
    ϕ[idx] = Float64(i)  # Linear: 1, 2, 3
end

# Set up a simple flow field (left to right)
for i in 1:N, j in 1:N
    idx = cartesian_to_linear_indices(i, j, mesh)
    u[idx] = 1.0  # Uniform flow in x-direction
    v[idx] = 0.0  # No flow in y-direction
end

println("Cell values ϕ:")
for j in N:-1:1
    for i in 1:N
        idx = cartesian_to_linear_indices(i, j, mesh)
        print("$(ϕ[idx])  ")
    end
    println()
end
println()

# Compute face values using different limiters
limiters = [VanLeer(), Minmod(), Superbee()]
limiter_names = ["VanLeer", "Minmod", "Superbee"]

for (limiter, name) in zip(limiters, limiter_names)
    println("Face values with $name limiter:")
    ϕ_f = compute_face_values(ϕ, mesh, u, v, limiter, 4.0)  # Top boundary value = 4

    for face_idx in 1:total_faces
        i1, j1, i2, j2, orientation = face_neighbors(mesh, face_idx)
        println("  Face $face_idx ($orientation, cells ($i1,$j1)-($i2,$j2)): $(round(ϕ_f[face_idx], digits=3))")
    end
    println()
end

println("=== Summary ===")
println("The face indexing system provides:")
println("1. Unique indices for all internal faces")
println("2. Consistent mapping between cells and faces")
println("3. Support for TVD schemes with various limiters")
println("4. Proper handling of boundary conditions")
println()
println("This enables efficient computation of convection operators")
println("for the lid-driven cavity solver using TVD schemes.")