using Test
using MEK9470
using LinearSolve
using LinearAlgebra
using Symbolics
using CairoMakie
# Analytical solution and RHS
@variables x, y
Dxx = Differential(x)^2
Dyy = Differential(y)^2
## Parameters
N = 32
mesh = CartesianMesh(N)
n = 2
m = 3
π = pi

# Grid
cellcenters = cell_center.(1:N^2, Ref(mesh))
xs, ys = first.(cellcenters), last.(cellcenters)

##
f(x, y) = sin(n*π*x) * sin(m*π*y)
L_f = Dxx(f(x, y)) + Dyy(f(x, y))
L_f = expand_derivatives(L_f)
L_f = build_function(L_f, x, y, expression=Val{false})

rhs = L_f.(xs, ys)
u = zeros(N^2)
r, M, b = laplacian(u, mesh, 1.0)

u_sol = M \ (rhs - b)
r, _ = laplacian(u_sol, mesh, 1.0)
r

u_exact = f.(xs, ys)
##
fig = Figure(size=(600, 1200))
ax = Axis(fig[1, 1], title="Solved")
ax_exact = Axis(fig[2, 1], title="Exact")
ax_error = Axis(fig[3, 1], title="Error")
ax_rhs = Axis(fig[4, 1], title="RHS")
hm = heatmap!(ax, xs, ys, u_sol)
hm_exact = heatmap!(ax_exact, xs, ys, u_exact)
hm_error = heatmap!(ax_error, xs, ys, abs.(u_sol - u_exact))
hm_rhs = heatmap!(ax_rhs, xs, ys, rhs)
Colorbar(fig[1, 2], hm)
Colorbar(fig[2, 2], hm_exact)
Colorbar(fig[3, 2], hm_error)
Colorbar(fig[4, 2], hm_rhs)
fig