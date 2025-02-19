using MEK9470
using LinearAlgebra
using CairoMakie
##Parameters
A = 0.0
B = 1.0
ρ = 1
Γ = 0.1
ϕ_A = 1.0
ϕ_B = 0.0
##
function assemble_system(discretization::AbstractDiscretization, F, D_w, D_e, N)
    @info "Discretization not implemented"
end

function assemble_system(discretization::CentralDiff, F, D_w, D_e, N)
    α_w = [0.0;ones(N-1)]
    α_e = [ones(N-1);0.0]
    a_w = (F/2 * (2 .- α_w) .+ D_w)
    a_e = (-F/2 * (2 .- α_e) .+ D_e)
    a_p = a_w + a_e
    M = Tridiagonal(-a_w[2:end], a_p, -a_e[1:end-1])
    b = zeros(N)
    b[1] = a_w[1]*ϕ_A
    b[end] = a_e[end]*ϕ_B
    return M, b
end

function assemble_system(discretization::Upwind, F, D_w, D_e, N)
    a_w = (max(0.0, F) .+ D_w)
    a_e = (max(-F, 0.0) .+ D_e)
    a_p = a_w + a_e
    M = Tridiagonal(-a_w[2:end], a_p, -a_e[1:end-1])
    b = zeros(N)
    b[1] = a_w[1]*ϕ_A
    b[end] = a_e[end]*ϕ_B
    return M, b
end

function assemble_system(discretization::Hybrid, F, D_w, D_e, N)
    α_w = [0.0;ones(N-1)]
    α_e = [ones(N-1);0.0]
    F_w = F/2 * (2 .- α_w)
    F_e = -F/2 * (2 .- α_e)
    a_w = max.(0.0, F_w .+ D_w, F)
    a_e = max.(-F, F_e .+ D_e, 0.0)
    a_p = a_w + a_e

    M = Tridiagonal(-a_w[2:end], a_p, -a_e[1:end-1])
    b = zeros(N)
    b[1] = a_w[1]*ϕ_A
    b[end] = a_e[end]*ϕ_B
    return M, b
end

function assemble_system(discretization::QUICK, F, D_w, D_e, N)
    α_w = 1.0
    α_e = 1.0

    D = D_w[2]

    a_w = D + 6/8 * α_w*F + 1/8*α_e*F + 3/8*(1-α_w)*F
    a_e = D - 3/8 * α_e*F - 6/8*(1-α_e)*F - 1/8*(1-α_w)*F
    a_ww = -1/8 * α_w*F
    a_p = a_w + a_e + a_ww

    M = zeros(N, N)
    b = zeros(N)
    for i in 3:N-1
        M[i, i-2] = -a_ww
        M[i, i-1] = -a_w
        M[i, i] = a_p
        M[i, i+1] = -a_e
    end

    #node 1
    M[1, 1+1] = -(D + 1/3*D - 3/8*F)
    M[1,1] = -M[1,1+1] + 8/3 *D + 1/4 *F + F
    b[1] = (8/3*D + 1/4*F + F)*ϕ_A

    #node 2
    M[2, 2-1] = -(D + 7/8*F + 1/8*F)
    M[2, 2+1] = -(D - 3/8*F)
    M[2,2] = -M[2, 2-1] - M[2, 2+1] - F/4
    b[2] = -1/4*F*ϕ_A

    #last node
    M[N, N-2] = 1/8*F
    M[N, N-1] = -(D + 1/3*D + 6/8*F)
    M[N, N] = -M[N, N-1] - M[N, N-2] + 8/3*D - F
    b[N] = (8/3*D - F)*ϕ_B

    return M, b
end

function solve_5_1(u=0.1, N=5, discretization::AbstractDiscretization = CentralDiff())
    F = ρ*u

    x, dx = create_1D_grid(A, B, N)

    δx = [dx/2; dx*ones(N-1); dx/2]
    D_w = Γ ./ δx[1:N]
    D_e = Γ ./ δx[end-N+1:end]


    M, b = assemble_system(discretization, F, D_w, D_e, N)
    @debug "M = $M"
    @debug "b = $b"

    ϕ = M\b
    return x, ϕ
end

function get_exact_solution(u)
    F = ρ*u
    x_fine = collect(LinRange(A, B, 100))
    ϕ_exact = (ϕ_B - ϕ_A) *(exp.(F.*x_fine./Γ) .- 1) ./ (exp.(F*B./Γ) .- 1) .+ ϕ_A
    return x_fine, ϕ_exact
end
##
using Logging
with_logger(ConsoleLogger(stderr, Logging.Debug)) do
    x, ϕ = solve_5_1(0.2, 5, QUICK())
    ϕ_exact, x_fine = get_exact_solution(0.2)
    fig = Figure()
    ax = Axis(fig[1,1])
    scatterlines!(ax, x, ϕ)
    lines!(ax, x_fine, ϕ_exact, color=:black, linestyle=:dash, label="Exact")
    display(fig)
    @show ϕ
end
##
# Create plots for each case
cases = ["i)", "ii)", "iii)", "bonus"]
us = [0.1, 2.5, 2.5, 0.2]
Ns = [5, 5, 20, 5]

fig = Figure(size=(800, 600))
axes = [Axis(fig[i÷3+1, mod1(i, 2)], title = "Case $(cases[i]): u = $(us[i]), N = $(Ns[i])") for i in 1:4]

discretizations = [CentralDiff(), Upwind(), Hybrid(), QUICK()]
discretization_names = ["Central", "Upwind", "Hybrid", "QUICK"]
colors = [:blue, :red, :green, :purple]

for (i, u) in enumerate(us)
    x_fine, ϕ_exact = get_exact_solution(u)
    lines!(axes[i], x_fine, ϕ_exact, color=:black, linestyle=:dash, label="Exact")
    
    for (j, (disc, name)) in enumerate(zip(discretizations, discretization_names))
        x, ϕ = solve_5_1(u, Ns[i], disc)
        scatterlines!(axes[i], x, ϕ, color=colors[j], label=name)
    end
end

# Add a single legend for all plots
fig[end+1,:] = Legend(fig, axes[1], "Methods", framevisible=false, merge=true, orientation=:horizontal)
Label(fig[0,:], "Example 5.1")
# Adjust layout
for ax in axes
    ax.xlabel = "x"
    ax.ylabel = "ϕ"
end

fig

##
save("plots/example_5-1_comparison.png", fig)