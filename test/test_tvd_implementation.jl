@testitem "Simple 1D Flow - Uniform Velocity" begin
    # Test case: 3x3 mesh, uniform flow in x-direction
    mesh = CartesianMesh(3)
    N = 3

    # Setup uniform flow in x-direction
    u = ones(N^2)
    v = zeros(N^2)

    # Linear variation in x: ϕ(i,j) = i
    ϕ = zeros(N^2)
    for i in 1:N, j in 1:N
        ϕ[cartesian_to_linear_indices(i, j, mesh)] = Float64(i)
    end

    # Test with different limiters
    limiters = [VanLeer(), Minmod(), Superbee()]

    for limiter in limiters
        bc = DomainBoundaryConditions(
            north=DirichletBC(0.0), south=DirichletBC(0.0),
            east=DirichletBC(0.0), west=DirichletBC(0.0)
        )
        ϕ_f = compute_face_values(ϕ, mesh, u, v, limiter, bc)

        # Face between cells (1,1) and (2,1): r = 0
        face_idx = get_face_index(mesh, 1, 1, :east)

        if limiter isa VanLeer
            # ψ(0) = 0 for VanLeer
            expected = 1.0 + 0.5 * 0.0 * 1.0
            @test ϕ_f[face_idx] ≈ expected atol = 1e-10
        end

        # Face between cells (2,1) and (3,1): r = 1
        face_idx = get_face_index(mesh, 2, 1, :east)

        if limiter isa VanLeer
            # ψ(1) = 1 for VanLeer
            expected = 2.0 + 0.5 * 1.0 * 1.0
            @test ϕ_f[face_idx] ≈ expected atol = 1e-10
        elseif limiter isa Minmod
            # ψ(1) = 1 for Minmod
            expected = 2.0 + 0.5 * 1.0 * 1.0
            @test ϕ_f[face_idx] ≈ expected atol = 1e-10
        end
    end
end

@testitem "Zero Velocity - Bounded Values" begin
    mesh = CartesianMesh(3)
    N = 3

    u = zeros(N^2)
    v = zeros(N^2)
    ϕ = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]

    limiter = VanLeer()
    bc = DomainBoundaryConditions(
        north=DirichletBC(0.0), south=DirichletBC(0.0),
        east=DirichletBC(0.0), west=DirichletBC(0.0)
    )
    ϕ_f = compute_face_values(ϕ, mesh, u, v, limiter, bc)

    @test all(isfinite.(ϕ_f))
    @test all(ϕ_f .>= minimum(ϕ))
    @test all(ϕ_f .<= maximum(ϕ))
end

@testitem "Single Cell Flow - Analytical Case" begin
    mesh = CartesianMesh(3)
    N = 3

    u = zeros(N^2)
    v = zeros(N^2)

    # Only cell (2,2) has velocity
    center_cell = cartesian_to_linear_indices(2, 2, mesh)
    u[center_cell] = 1.0

    # Step function in ϕ
    ϕ = zeros(N^2)
    ϕ[cartesian_to_linear_indices(1, 2, mesh)] = 1.0
    ϕ[cartesian_to_linear_indices(2, 2, mesh)] = 2.0
    ϕ[cartesian_to_linear_indices(3, 2, mesh)] = 3.0

    limiter = VanLeer()
    bc = DomainBoundaryConditions(
        north=DirichletBC(0.0), south=DirichletBC(0.0),
        east=DirichletBC(0.0), west=DirichletBC(0.0)
    )
    ϕ_f = compute_face_values(ϕ, mesh, u, v, limiter, bc)

    # Test the east face of cell (2,2)
    # Expected: ϕ_face = 2 + 0.5 * 1 * 1 = 2.5
    face_idx = get_face_index(mesh, 2, 2, :east)
    expected = 2.0 + 0.5 * 1.0 * 1.0
    @test ϕ_f[face_idx] ≈ expected atol = 1e-10
end

@testitem "Reverse Flow Test" begin
    mesh = CartesianMesh(3)
    N = 3

    u = -ones(N^2)  # Negative flow
    v = zeros(N^2)

    # Linear variation: ϕ(i,j) = i
    ϕ = zeros(N^2)
    for i in 1:N, j in 1:N
        ϕ[cartesian_to_linear_indices(i, j, mesh)] = Float64(i)
    end

    limiter = VanLeer()
    bc = DomainBoundaryConditions(
        north=DirichletBC(0.0), south=DirichletBC(0.0),
        east=DirichletBC(0.0), west=DirichletBC(0.0)
    )
    ϕ_f = compute_face_values(ϕ, mesh, u, v, limiter, bc)

    # For negative flow: ϕ_face = 2 + 0.5 * 1 * (-1) = 1.5
    face_idx = get_face_index(mesh, 1, 1, :east)
    expected = 2.0 + 0.5 * 1.0 * (-1.0)
    @test ϕ_f[face_idx] ≈ expected atol = 1e-10
end

@testitem "Boundary Conditions Test" begin
    mesh = CartesianMesh(3)
    N = 3

    u = zeros(N^2)
    v = ones(N^2)  # Upward flow
    ϕ = 2.0 * ones(N^2)

    # Test Dirichlet BC
    bc_dirichlet = DomainBoundaryConditions(
        north=DirichletBC(5.0), south=DirichletBC(0.0),
        east=DirichletBC(0.0), west=DirichletBC(0.0)
    )
    ϕ_f_d = compute_face_values(ϕ, mesh, u, v, VanLeer(), bc_dirichlet)
    @test all(isfinite.(ϕ_f_d))

    # Test Neumann BC
    bc_neumann = DomainBoundaryConditions(
        north=NeumannBC(0.0), south=NeumannBC(0.0),
        east=NeumannBC(0.0), west=NeumannBC(0.0)
    )
    ϕ_f_n = compute_face_values(ϕ, mesh, u, v, VanLeer(), bc_neumann)
    @test all(isfinite.(ϕ_f_n))

    # Results should be different for different BCs
    @test ϕ_f_d != ϕ_f_n
end

@testitem "Limiter Function Comparison" begin
    mesh = CartesianMesh(3)
    N = 3

    u = ones(N^2)
    v = zeros(N^2)
    ϕ = [1.0, 3.0, 4.0, 1.0, 3.0, 4.0, 1.0, 3.0, 4.0]

    limiters = [VanLeer(), Minmod(), Superbee()]
    bc = DomainBoundaryConditions(
        north=DirichletBC(0.0), south=DirichletBC(0.0),
        east=DirichletBC(0.0), west=DirichletBC(0.0)
    )
    results = [compute_face_values(ϕ, mesh, u, v, limiter, bc) for limiter in limiters]

    # Results should be different for different limiters
    @test results[1] != results[2]
    @test results[1] != results[3]
    @test results[2] != results[3]

    # All should be finite and bounded
    for result in results
        @test all(isfinite.(result))
        @test all(result .>= minimum(ϕ) - 1e-10)
        @test all(result .<= maximum(ϕ) + 1e-10)
    end
end

@testitem "r-ratio Calculation" begin
    mesh = CartesianMesh(3)
    ϕ = [1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0]
    bc = DomainBoundaryConditions(
        north=DirichletBC(0.0), south=DirichletBC(0.0),
        east=DirichletBC(0.0), west=DirichletBC(0.0)
    )

    # Test r calculation for flow from (1,1) to (2,1)
    r = compute_r_ratio(ϕ, mesh, 1, 1, 2, 1, Vertical(), bc)
    @test r ≈ 0.0 atol = 1e-10

    # Test r calculation for flow from (2,1) to (3,1)
    r = compute_r_ratio(ϕ, mesh, 2, 1, 3, 1, Vertical(), bc)
    @test r ≈ 1.0 atol = 1e-10
end

@testitem "Multiple Dispatch Functionality" begin
    # Test that our dispatch system works correctly
    @test face_velocity(1.0, 2.0, Vertical()) == 1.0
    @test face_velocity(1.0, 2.0, Horizontal()) == 2.0

    @test flow_direction(1.0) isa Positive
    @test flow_direction(-1.0) isa Negative
    @test flow_direction(0.0) isa Positive

    @test upwind_cell(1, 2, Positive()) == 1
    @test upwind_cell(1, 2, Negative()) == 2
    @test downwind_cell(1, 2, Positive()) == 2
    @test downwind_cell(1, 2, Negative()) == 1
end

@testitem "Convection Operators" begin
    mesh = CartesianMesh(3)
    N = 3

    u = ones(N^2)
    v = zeros(N^2)
    ϕ = ones(N^2)

    limiter = VanLeer()
    bc = DomainBoundaryConditions(
        north=DirichletBC(0.0), south=DirichletBC(0.0),
        east=DirichletBC(0.0), west=DirichletBC(0.0)
    )

    M_u, M_v = get_convection_operators(mesh, ϕ, u, v, limiter, bc)

    @test size(M_u) == (N^2, N^2)
    @test size(M_v) == (N^2, N^2)
    @test all(isfinite.(M_u))
    @test all(isfinite.(M_v))
end