@testitem "Total faces calculation" begin
    # Test with a small 3x3 mesh
    mesh = CartesianMesh(3)
    N = 3

    total_faces = get_total_faces(mesh)
    expected = 2 * N * (N - 1)  # 2 * 3 * 2 = 12
    @test total_faces == expected
end

@testitem "Face indexing consistency" begin
    mesh = CartesianMesh(3)
    N = 3

    # Test that east/west faces are consistent
    for i in 1:N-1, j in 1:N
        east_idx = get_face_index(mesh, i, j, :east)
        west_idx = get_face_index(mesh, i + 1, j, :west)
        @test east_idx == west_idx
    end

    # Test that north/south faces are consistent
    for i in 1:N, j in 1:N-1
        north_idx = get_face_index(mesh, i, j, :north)
        south_idx = get_face_index(mesh, i, j + 1, :south)
        @test north_idx == south_idx
    end
end

@testitem "Face neighbors" begin
    mesh = CartesianMesh(3)
    N = 3

    # Test a few specific faces
    i1, j1, i2, j2, orientation = face_neighbors(mesh, 1)
    @test orientation == :vertical
    @test abs(i2 - i1) == 1 && j1 == j2

    # Test horizontal face
    vertical_faces = N * (N - 1)
    i1, j1, i2, j2, orientation = face_neighbors(mesh, vertical_faces + 1)
    @test orientation == :horizontal
    @test i1 == i2 && abs(j2 - j1) == 1
end

@testitem "Boundary assertions" begin
    mesh = CartesianMesh(3)
    N = 3

    # Test that boundary faces throw errors appropriately
    @test_throws AssertionError get_face_index(mesh, N, 1, :east)  # No east face for rightmost
    @test_throws AssertionError get_face_index(mesh, 1, 1, :west)  # No west face for leftmost
    @test_throws AssertionError get_face_index(mesh, 1, N, :north) # No north face for topmost
    @test_throws AssertionError get_face_index(mesh, 1, 1, :south) # No south face for bottommost
end

@testitem "Uniform field test" begin
    mesh = CartesianMesh(4)
    N = 4

    # Create simple test data
    ϕ = ones(N^2)  # Uniform field
    u = zeros(N^2)  # No flow
    v = zeros(N^2)
    limiter = VanLeer()

    ϕ_f = compute_face_values(ϕ, mesh, u, v, limiter, DomainBoundaryConditions())
    # For uniform field, all face values should be 1.0
    @test all(ϕ_f .≈ 1.0)
end

@testitem "Linear field test" begin
    mesh = CartesianMesh(4)
    N = 4

    # Create simple test data
    ϕ = ones(N^2)  # Start with uniform field
    u = zeros(N^2)  # No flow
    v = zeros(N^2)
    limiter = VanLeer()

    # Create a linear field in x-direction
    for i in 1:N, j in 1:N
        idx = cartesian_to_linear_indices(i, j, mesh)
        ϕ[idx] = i  # Linear in x
    end

    ϕ_f = compute_face_values(ϕ, mesh, u, v, limiter, DomainBoundaryConditions())

    # For linear field with no flow, face values should be averages
    # This is a basic sanity check
    @test length(ϕ_f) == get_total_faces(mesh)
    @test all(isfinite.(ϕ_f))
end