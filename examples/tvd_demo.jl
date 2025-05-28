using MEK9470

println("=== TVD Scheme Demonstration ===")

# Create test case with step function
mesh = CartesianMesh(5)
N = 5

u = ones(N^2)  # Uniform flow in x-direction
v = zeros(N^2)

# Step function: ϕ = 1 for i ≤ 2, ϕ = 3 for i > 2
ϕ = [i ≤ 2 ? 1.0 : 3.0 for i in 1:N for j in 1:N]

println("Initial ϕ field (step function):")
for j in N:-1:1, i in 1:N
    print("$(ϕ[cartesian_to_linear_indices(i, j, mesh)]) ")
    i == N && println()
end
println()

# Test different limiters with Dirichlet BC
limiters = [
    ("Upwind (UD)", UD()),
    ("Van Leer", VanLeer()),
    ("Minmod", Minmod()),
    ("Superbee", Superbee()),
    ("UMIST", UMIST())
]

bc = DirichletBC(0.0)

for (name, limiter) in limiters
    println("=== $name Limiter ===")

    ϕ_f = compute_face_values(ϕ, mesh, u, v, limiter, bc)

    # Show face values at the critical interface (i=2 to i=3)
    println("Face values at step interface:")
    for j in 1:N
        face_idx = get_face_index(mesh, 2, j, :east)
        println("  Face (2,$j)→(3,$j): ϕ_f = $(round(ϕ_f[face_idx], digits=3))")
    end

    # Show r-ratios and limiter values
    println("r-ratios at interface:")
    for j in 1:N
        r = compute_r_ratio(ϕ, mesh, 2, j, 3, j, Vertical(), bc)
        ψ_r = limiter(r)
        println("  Face (2,$j)→(3,$j): r = $(round(r, digits=3)), ψ(r) = $(round(ψ_r, digits=3))")
    end
    println()
end

# Demonstrate limiter function behavior
println("=== Limiter Function Behavior ===")
r_values = [0.0, 0.5, 1.0, 1.5, 2.0, 5.0]

println("r\t\tUD\tVanLeer\tMinmod\tSuperbee\tUMIST")
println("-"^60)

for r in r_values
    vals = [round(limiter(r), digits=3) for (_, limiter) in limiters]
    println("$r\t\t$(join(vals, "\t"))")
end

println("\n=== TVD Properties ===")
properties = [
    "UD (Upwind): ψ(r) = 0 → First-order upwind (most diffusive)",
    "Van Leer: ψ(r) = (r + |r|)/(1 + |r|) → Smooth, monotonic",
    "Minmod: ψ(r) = max(0, min(1, r)) → Most compressive",
    "Superbee: ψ(r) = max(0, min(1, 2r), min(2, r)) → Least diffusive",
    "UMIST: ψ(r) = max(0, min(2r, (1+3r)/4, (3+r)/4, 2)) → Good compromise"
]

for (i, prop) in enumerate(properties)
    println("$i. $prop")
end

println("\nTVD constraint: 0 ≤ ψ(r) ≤ min(2r, 2) for r > 0")

# Demonstrate boundary condition effects
println("\n=== Boundary Condition Comparison ===")
bc_dirichlet = DirichletBC(5.0)
bc_neumann = NeumannBC(0.0)

ϕ_f_d = compute_face_values(ϕ, mesh, u, v, VanLeer(), bc_dirichlet)
ϕ_f_n = compute_face_values(ϕ, mesh, u, v, VanLeer(), bc_neumann)

println("Van Leer with different boundary conditions:")
println("Dirichlet BC (top = 5.0): max face value = $(round(maximum(ϕ_f_d), digits=3))")
println("Neumann BC (∇ϕ = 0.0): max face value = $(round(maximum(ϕ_f_n), digits=3))")