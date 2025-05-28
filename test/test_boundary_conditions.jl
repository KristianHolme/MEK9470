@testitem "Domain Boundary Conditions Construction" begin
    # Test basic construction
    bc = DomainBoundaryConditions(
        north=DirichletBC(1.0),
        south=DirichletBC(0.0),
        east=NeumannBC(0.0),
        west=DirichletBC(0.5)
    )

    @test bc.north isa DirichletBC
    @test bc.south isa DirichletBC
    @test bc.east isa NeumannBC
    @test bc.west isa DirichletBC

    @test bc.north.value == 1.0
    @test bc.south.value == 0.0
    @test bc.east.gradient == 0.0
    @test bc.west.value == 0.5
end

@testitem "Lid-Driven Cavity Convenience Functions" begin
    # Test u-velocity BC
    u_bc = lid_driven_cavity_u_bc(2.0)
    @test u_bc.north isa DirichletBC
    @test u_bc.north.value == 2.0
    @test u_bc.south.value == 0.0
    @test u_bc.east.value == 0.0
    @test u_bc.west.value == 0.0

    # Test v-velocity BC
    v_bc = lid_driven_cavity_v_bc()
    @test all(bc -> bc isa DirichletBC && bc.value == 0.0,
        [v_bc.north, v_bc.south, v_bc.east, v_bc.west])
end

@testitem "Boundary Side Detection" begin
    N = 4

    # Test corner and edge cases
    @test get_boundary_side(1, 1, N) == South()  # Bottom-left corner (south takes precedence)
    @test get_boundary_side(N, 1, N) == South()  # Bottom-right corner
    @test get_boundary_side(1, N, N) == North()  # Top-left corner
    @test get_boundary_side(N, N, N) == North()  # Top-right corner

    # Test edges
    @test get_boundary_side(2, 1, N) == South()  # Bottom edge
    @test get_boundary_side(2, N, N) == North()  # Top edge
    @test get_boundary_side(1, 2, N) == West()   # Left edge
    @test get_boundary_side(N, 2, N) == East()   # Right edge

    # Test interior
    @test get_boundary_side(2, 2, N) === nothing  # Interior cell
    @test get_boundary_side(3, 3, N) === nothing  # Interior cell
end

@testitem "Boundary Condition Dispatch" begin
    bc = DomainBoundaryConditions(
        north=DirichletBC(1.0),
        south=DirichletBC(0.0),
        east=NeumannBC(0.5),
        west=DirichletBC(-1.0)
    )

    @test get_bc(bc, North()) == bc.north
    @test get_bc(bc, South()) == bc.south
    @test get_bc(bc, East()) == bc.east
    @test get_bc(bc, West()) == bc.west

    @test get_bc(bc, North()).value == 1.0
    @test get_bc(bc, South()).value == 0.0
    @test get_bc(bc, East()).gradient == 0.5
    @test get_bc(bc, West()).value == -1.0
end

@testitem "TVD with Domain Boundary Conditions" begin
    N = 3
    mesh = CartesianMesh(N)

    # Create simple test fields
    u = ones(N^2) * 0.1
    v = zeros(N^2)
    ϕ = collect(1.0:N^2)

    # Test with lid-driven cavity BCs
    u_bc = lid_driven_cavity_u_bc(1.0)
    limiter = VanLeer()

    # Should not throw errors
    ϕ_f = compute_face_values(ϕ, mesh, u, v, limiter, u_bc)
    @test length(ϕ_f) == get_total_faces(mesh)
    @test all(isfinite, ϕ_f)

    # Test convection operators
    M_u, M_v = get_convection_operators(mesh, ϕ, u, v, limiter, u_bc)
    @test size(M_u) == (N^2, N^2)
    @test size(M_v) == (N^2, N^2)
end

@testitem "Different BCs on Different Sides" begin
    N = 3
    mesh = CartesianMesh(N)

    u = ones(N^2) * 0.1
    v = zeros(N^2)
    ϕ = ones(N^2)  # Uniform field
    limiter = UD()  # Upwind for simplicity

    # Create BCs with different values on different sides
    bc1 = DomainBoundaryConditions(
        north=DirichletBC(1.0), south=DirichletBC(0.0),
        east=DirichletBC(0.0), west=DirichletBC(0.0)
    )

    bc2 = DomainBoundaryConditions(
        north=DirichletBC(0.0), south=DirichletBC(0.0),
        east=DirichletBC(0.0), west=DirichletBC(0.0)
    )

    ϕ_f1 = compute_face_values(ϕ, mesh, u, v, limiter, bc1)
    ϕ_f2 = compute_face_values(ϕ, mesh, u, v, limiter, bc2)

    # Results should be different due to different boundary conditions
    @test !(ϕ_f1 ≈ ϕ_f2)

    # But both should be finite and reasonable
    @test all(isfinite, ϕ_f1)
    @test all(isfinite, ϕ_f2)
end