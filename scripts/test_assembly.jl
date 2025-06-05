using DrWatson
@quickactivate :MEK9470
using LinearAlgebra
using LinearSolve
using SparseArrays
using CairoMakie
using Dates
##
mesh = CartesianMesh(32)

conv = ConvectionOperator(mesh)
uvp = 1:mesh.N^2*3 |> collect .|> Float64
update_convection_operators!(conv, mesh, uvp, VanLeer(), lid_driven_cavity_u_bc(1.0), lid_driven_cavity_v_bc())

ψ = VanLeer()
ldcprob = LDCProblem(mesh, 1.0, ψ, 1.0)

M
sparse(b)
##
M, b = construct_linear_system(ldcprob, uvp)
linprob = LinearProblem(M, b)
sol = solve(linprob)
uvp = sol.u
lid_driven_cavity_plot(sol.u, mesh, title="Linear system solution")
##
using NonlinearSolve
mesh = CartesianMesh(32)
ψ = VanLeer()
ldcprob = LDCProblem(mesh, 0.001, ψ, 1.0)
function residuals(uvp, ldcprob)
    M, b = construct_linear_system(ldcprob, uvp)
    return M * uvp - b
end

# Test nonlinear solve with AD-compatible types
uvp = zeros(mesh.N^2 * 3) .+ 0.5
nlprob = NonlinearProblem(residuals, uvp, ldcprob)
sol = solve(nlprob)
fig = lid_driven_cavity_plot(sol.u, mesh, title="Converged Nonlinear Solution")
save(plotsdir("64x64_solution_$(Dates.format(now(), "yyyy-mm-dd_HH-MM")).png"), fig)
norm(residuals(sol.u, ldcprob))
M, b = construct_linear_system(ldcprob, sol.u)
norm(M)
norm(b)
# Plot the converged solution
plot_velocity_profiles(sol.u, mesh)
plot_streamfunction(sol.u, mesh)
