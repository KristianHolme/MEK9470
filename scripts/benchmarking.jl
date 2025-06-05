using DrWatson
@quickactivate :MEK9470
using LinearAlgebra
using LinearSolve
using SparseArrays
using CairoMakie
using Dates
using NonlinearSolve
using BenchmarkTools
using JLD2
##
mesh = CartesianMesh(16)
ψ = VanLeer()
ldcprob = LDCProblem(mesh, 0.1, ψ, 1.0)
uvp = zeros(mesh.N^2 * 3) .+ 0.5
M, b = construct_linear_system(ldcprob, uvp)
M
b


# Test nonlinear solve with AD-compatible types
uvp = zeros(mesh.N^2 * 3) .+ 0.5
sol = solve(nlprob, NewtonRaphson())
##
using SparseConnectivityTracer, ADTypes
using NonlinearSolve.NonlinearSolveBase: RelNormSafeBestTerminationMode
uvp = zeros(mesh.N^2 * 3) .+ 0.5
##
mesh = CartesianMesh(64)
ψ = QUICKlimiter()
##
Re = 1
ν = 1 / Re
ldcprob = LDCProblem(mesh, ν, ψ)
jac_sparsity = ADTypes.jacobian_sparsity(u -> residuals(u, ldcprob), uvp, TracerLocalSparsityDetector())

ff = NonlinearFunction(residuals; jac_prototype=jac_sparsity)
termination_condition = RelNormSafeBestTerminationMode(norm)
tol = 1e-4
nlprob_sparse = NonlinearProblem(ff, uvp, ldcprob; abstol=tol, reltol=tol, termination_condition=termination_condition)
##
@profview sol = solve(nlprob_sparse, NewtonRaphson())
if sol.retcode == ReturnCode.Success
    uvp = sol.u
    dir = datadir("solutions", "N=$(mesh.N)", "Re=$(Re)")
    mkpath(dir)
    jldsave(joinpath(dir, "$(string(ψ)).jld2"); uvp)
end
sol.retcode == ReturnCode.Success
##
fig = lid_driven_cavity_plot(sol.u, mesh, title="N=$(mesh.N), Re=$(Re), $(string(ψ))", xlims=(0, 1), ylims=(0, 1),
    show_streamlines=true)
safesave(plotsdir("ldc", "N=$(mesh.N)", "Re=$(Re)", "$(string(ψ)).png"), fig)
extrema(uvp)
norm(residuals(uvp, ldcprob))
##
##
@time solve(nlprob, NewtonRaphson());
@time solve(nlprob_sparse);
@time solve(nlprob_sparse, NewtonRaphson());
@time solve(nlprob_sparse, NewtonRaphson(linsolve=KrylovJL_GMRES()));
@time solve(nlprob, NewtonRaphson(linsolve=KLUFactorization()));
##