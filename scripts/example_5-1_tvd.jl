using MEK9470
using LinearAlgebra
using CairoMakie
using Logging
##Parameters
A = 0.0
B = 1.0
ρ = 1
Γ = 0.1
ϕ_A = 1.0
ϕ_B = 0.0
##

function r_e(ϕ, P)
    @assert P ∈ eachindex(ϕ)
    E = P+1
    W = P-1
    if P == firstindex(ϕ)
        r = 2*(ϕ[P] - ϕ_A) / (ϕ[E] - ϕ[P])
        r > 0 || @logmsg Logging.LogLevel(-1) "r_e = $r is negative for P = $P, ϕ = $ϕ"
        return r
    elseif P == lastindex(ϕ)
        r = (ϕ[P] - ϕ[W]) / (2*(ϕ_B - ϕ[P]))
        r > 0 || @logmsg Logging.LogLevel(-1) "r_e = $r is negative for P = $P, ϕ = $ϕ"
        return r
    else
        r = (ϕ[P] - ϕ[W]) / (ϕ[E] - ϕ[P])
        r > 0 || @logmsg Logging.LogLevel(-1) "r_e = $r is negative for P = $P, ϕ = $ϕ"
        return r
    end
end

function r_w(ϕ, P)
    @assert P ∈ eachindex(ϕ)
    W = P-1
    WW = P-2
    if P == firstindex(ϕ)
        @error "Dont know how to handle ϕ[WW] for P = $P"
    elseif P == firstindex(ϕ) + 1
        r = 2*(ϕ[W] - ϕ_A) / (ϕ[P] - ϕ[W])
        r > 0 || @logmsg Logging.LogLevel(-1) "r_w = $r is negative for P = $P, ϕ = $ϕ"
        return r
    else
        r = (ϕ[W] - ϕ[WW]) / (ϕ[P] - ϕ[W])
        r > 0 || @logmsg Logging.LogLevel(-1) "r_w = $r is negative for P = $P, ϕ = $ϕ"
        return r
    end
end

function update_b!(b, ϕ, D, F, Discretization)      
    ψ(r) = apply_limiter(Discretization.limiter, r)
    b[1] = (2D+F)*ϕ_A - F/2 * ψ(r_e(ϕ, firstindex(ϕ)))*(ϕ[2] - ϕ[1])
    for P in 2:lastindex(ϕ)-1
        E = P+1
        W = P-1
        b[P] = - F/2 * ψ(r_e(ϕ, P))*(ϕ[E] - ϕ[P]) 
               + F/2 * ψ(r_w(ϕ, P))*(ϕ[P] - ϕ[W])
    end
    b[end] = (-F + 2D)*ϕ_B + F/2 * ψ(r_w(ϕ, lastindex(ϕ)))*(ϕ[end] - ϕ[end-1]) 
    return b
end

function solve_tvd(N, u, Discretization)
    F = ρ*u
    x, dx = create_1D_grid(A, B, N)
    D = Γ / dx
    a_w = (D+F) .* ones(N-1)
    a_e = D .* ones(N-1)
    a_e[1] = D
    a_p = (2D+F) .* ones(N)
    a_p[1] = 2D + F + D
    a_p[end] = 2D + D
    M = Tridiagonal(-a_w, a_p, -a_e)
    @show M
    ϕ = ϕ_A .* (1 .- x) + ϕ_B .* x
    ϕ_prev = ϕ .* 0
    b = zeros(N)
    update_b!(b, ϕ, D, F, Discretization)
    @show b
    error = norm(ϕ - ϕ_prev)
    its = 0
    while error > 1e-6 && its < 1000
        its += 1
        @debug "Iteration $its"
        copyto!(ϕ_prev, ϕ)
        update_b!(b, ϕ, D, F, Discretization)
        ϕ .= M\b
        error = norm(ϕ - ϕ_prev)
    end
    return ϕ, x, its
end

function get_exact_solution(u)
    F = ρ*u
    x_fine = collect(LinRange(A, B, 100))
    ϕ_exact = (ϕ_B - ϕ_A) *(exp.(F.*x_fine./Γ) .- 1) ./ (exp.(F*B./Γ) .- 1) .+ ϕ_A
    return x_fine, ϕ_exact
end
##
Discretization = TVD(UD())
ϕ, x, its = solve_tvd(20, 2.5, Discretization)
its

##
# Create plots for each case
cases = ["i)", "i.2)", "ii)", "iii)"]
us = [0.1, 0.1,  2.5, 2.5]
Ns = [5, 20, 5, 20]

fig = Figure(size=(800, 600))
axes = [Axis(fig[i÷3+1, mod1(i, 2)], title = "Case $(cases[i]): u = $(us[i]), N = $(Ns[i])") for i in eachindex(cases)]

idxs = 1:7
limiters = [UD(),
            VanLeer(),
            VanAlbada(),
            Minmod(),
            Superbee(),
            Sweby(),
            UMIST()][idxs]
limiter_names = ["UD",
                 "VanLeer",
                 "VanAlbada",
                 "Minmod",
                 "Superbee",
                 "Sweby",
                 "UMIST"][idxs]

for (i, u) in enumerate(us)
    @debug "case $i"
    x_fine, ϕ_exact = get_exact_solution(u)
    lines!(axes[i], x_fine, ϕ_exact, color=:black, linestyle=:dash, label="Exact")
    
    for (j, (limiter, name)) in enumerate(zip(limiters, limiter_names))
        discretization = TVD(limiter)
        @debug "Solving for $name"
        ϕ, x, its = solve_tvd(Ns[i], u, discretization)
        scatterlines!(axes[i], x, ϕ, label=name)
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