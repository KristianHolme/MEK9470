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
function solve_5_1(u=0.1, N=5)
    F = ρ*u

    x, dx = create_1D_grid(A, B, N)
    α_w = [0.0;ones(N-1)]
    α_e = [ones(N-1);0.0]
    δx = [dx/2; dx*ones(N-1); dx/2]
    D_w = Γ ./ δx[1:N]
    D_e = Γ ./ δx[end-N+1:end]


    a_w = (F/2 * (2 .- α_w) .+ D_w)
    a_e = (-F/2 * (2 .- α_e) .+ D_e)
    a_p = a_w + a_e

    M = Tridiagonal(-a_w[2:end], a_p, -a_e[1:end-1])
    b = zeros(N)
    b[1] = a_w[1]*ϕ_A
    b[end] = a_e[end]*ϕ_B
    b

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

cases = ["i)", "ii)", "iii)"]
us = [0.1, 2.5, 2.5]
N = [5, 5, 20]

fig = Figure()
ax = Axis(fig[1,1], title = "Example 5.1")
for i in 1:3
    x, ϕ = solve_5_1(us[i], N[i])
    x_fine, ϕ_exact = get_exact_solution(us[i])

    scatterlines!(ax, x, ϕ, label = "u = $(us[i]), N = $(N[i])")
    lines!(ax, x_fine, ϕ_exact, linestyle = :dash, color = :red, label = "Exact, u = $(us[i])")
end
axislegend(ax, merge = true)
fig
save("plots/example_5-1.png", fig)