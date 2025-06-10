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
using SparseConnectivityTracer, ADTypes
using NonlinearSolve.NonlinearSolveBase: RelNormSafeBestTerminationMode
##
limiters = [Sweby(), VanLeer(), VanAlbada(), Minmod(), Superbee(), UMIST(), QUICKlimiter(), UD(), CD()]
##
N = 128
cmesh = CartesianMesh(N)
ψ = limiters[4]
uvp = zeros(cmesh.N^2 * 3) .+ 0.5

##
# 1, 10, 20, 40, 70, 100, 150, 200, 250, 350, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000
Re = 1500
ν = 1 / Re
ldcprob = LDCProblem(cmesh, ν, ψ)
jac_sparsity = ADTypes.jacobian_sparsity(u -> residuals(u, ldcprob), uvp, TracerLocalSparsityDetector())

ff = NonlinearFunction(residuals; jac_prototype=jac_sparsity)
tol = 1e-4
nlprob_sparse = NonlinearProblem(ff, uvp, ldcprob; abstol=tol, reltol=tol, termination_condition=RelNormSafeBestTerminationMode(norm))
##
# @profview sol = solve(nlprob_sparse, NewtonRaphson())
t_start = time()
sol = solve(nlprob_sparse, NewtonRaphson());
t_end = time()
t_elapsed = t_end - t_start
if sol.retcode == ReturnCode.Success
    uvp = sol.u
end;
sol.retcode
##
dir = datadir("solutions", "N=$(cmesh.N)", "Re=$(Re)")
isdir(dir) || mkpath(dir)
res_norm = norm(residuals(uvp, ldcprob))
jldsave(joinpath(dir, "$(string(ψ)).jld2"); uvp, res_norm, t_elapsed)

# Create and save plot
fig = lid_driven_cavity_plot(sol.u, cmesh, title="N=$(cmesh.N), Re=$(Re), $(string(ψ))",
    xlims=(0, 1), ylims=(0, 1), velocity_contours=true)
display(fig)
safesave(plotsdir("ldc", "N=$(cmesh.N)", "Re=$(Re)", "$(string(ψ)).png"), fig)
## For Re = 250, 1000
stream = compute_streamfunction(sol.u, cmesh)
fig = Figure(size=(600, 600))
ax = Axis(fig[1, 1], title="Streamfunction, Re = $(Re), ψ = $(string(ψ))", xlabel="x", ylabel="y", aspect=DataAspect())
stream_mat = reshape(stream, cmesh.N, cmesh.N)
coords = cell_center.(1:N, 1:N, Ref(cmesh))
xs = getindex.(coords, Ref(1))
ys = getindex.(coords, Ref(2))
hm = heatmap!(ax, xs, ys, stream_mat, colormap=:coolwarm)
I = argmin(stream_mat)
i, j = I[1], I[2]
ccx, ccy = cell_center(i, j, cmesh)
sc = scatter!(ax, [ccx], [ccy], color=:magenta, markersize=10, label="Minimum")
text!(ax, ccx + 0.01, ccy + 0.01, text=string("($(round(ccx, digits=3)), $(round(ccy, digits=3)))"), fontsize=12)
Colorbar(fig[1, 2], hm)

display(fig)
safesave(plotsdir("ldc", "N=$(cmesh.N)", "Re=$(Re)", "$(string(ψ))_streamfunction.png"), fig)







##
# @time solve(nlprob, NewtonRaphson());
# @time solve(nlprob_sparse);
# @time solve(nlprob_sparse, NewtonRaphson());
# @time solve(nlprob_sparse, NewtonRaphson(linsolve=KrylovJL_GMRES()));
# @time solve(nlprob, NewtonRaphson(linsolve=KLUFactorization()));
##