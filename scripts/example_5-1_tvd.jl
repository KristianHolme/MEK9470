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
ϕ_B = 0.3
##

function r_e(ϕ, P)
    @assert P ∈ eachindex(ϕ)
    E = P+1
    W = P-1
    if P == firstindex(ϕ)
        r = 2*(ϕ[P] - ϕ_A) / (ϕ[E] - ϕ[P])
    elseif P == lastindex(ϕ)
        r = (ϕ[P] - ϕ[W]) / (2*(ϕ_B - ϕ[P]))
    else
        r = (ϕ[P] - ϕ[W]) / (ϕ[E] - ϕ[P])
    end
    r > 0 || @logmsg Logging.LogLevel(-1) "r_e = $r is negative for P = $P"
    return r
end

function r_w(ϕ, P)
    @assert P ∈ eachindex(ϕ)
    W = P-1
    WW = P-2
    if P == firstindex(ϕ)
        @error "Dont know how to handle ϕ[WW] for P = $P"
        r = nothing
    elseif P == firstindex(ϕ) + 1
        r = 2*(ϕ[W] - ϕ_A) / (ϕ[P] - ϕ[W])
    else
        r = (ϕ[W] - ϕ[WW]) / (ϕ[P] - ϕ[W])
    end
    r > 0 || @logmsg Logging.LogLevel(-1) "r_w = $r is negative for P = $P"
    return r
end

function update_b!(b, ϕ, D, F, Discretization)      
    ψ(r) = apply_limiter(Discretization.limiter, r)
    b[1] = (D+F)*ϕ_A - F/2 * ψ(r_e(ϕ, firstindex(ϕ)))*(ϕ[2] - ϕ[1])
    for P in 2:(lastindex(ϕ)-1)
        E = P+1
        W = P-1
        b[P] = - F/2 * ψ(r_e(ϕ, P))*(ϕ[E] - ϕ[P]) 
               + F/2 * ψ(r_w(ϕ, P))*(ϕ[P] - ϕ[W])
    end
    b[end] = (D-F)*ϕ_B + F/2 * ψ(r_w(ϕ, lastindex(ϕ)))*(ϕ[end] - ϕ[end-1]) 
    return b
end

function solve_tvd(N, u, Discretization)
    F = ρ*u
    x, dx = create_1D_grid(A, B, N)
    D = Γ / dx
    a_w = (D+F) .* ones(N-1)
    a_e = D .* ones(N-1)
    a_p = (2D+F) .* ones(N)
    a_p[end] = 2D
    M = Tridiagonal(-a_w, a_p, -a_e)
    @debug "M = $M"
    ϕ = ϕ_A .* (1 .- x) + ϕ_B .* x
    ϕ_prev = zeros(N)
    b = zeros(N)
    update_b!(b, ϕ, D, F, Discretization)
    @debug "b = $b"
    error = norm(ϕ - ϕ_prev)
    its = 0
    while error > 1e-10 && its < 1000
        its += 1
        @debug "Iteration $its"
        copyto!(ϕ_prev, ϕ)
        update_b!(b, ϕ, D, F, Discretization)
        ϕ .= M\b
        error = norm(ϕ - ϕ_prev)
    end
    return ϕ, x, its
end

function get_exact_solution(u;x=collect(LinRange(A, B, 100)))
    F = ρ*u
    L = B-A
    ϕ_exact = (ϕ_B - ϕ_A) * (exp.(F.*x./Γ) .- 1) ./ (exp.(F*L./Γ) .- 1) .+ ϕ_A
    return ϕ_exact, x
end

function error(ϕ, ϕ_exact, dx)
    # L2 error is the square root of the integral of the squared difference
    integrand = (ϕ - ϕ_exact).^2
    error = 0.0
    for i in 1:length(ϕ)-1
        error += dx * (integrand[i] + integrand[i+1]) / 2
    end
    return sqrt(error)  # Need to take the square root for L2 norm
end
##
using Logging
Discretization = TVD(UMIST())
# Discretization = TVD(UD())
u = 0.1
N = 8
with_logger(ConsoleLogger(stderr, Logging.Info)) do
    ϕ, x, its = solve_tvd(N, u, Discretization);
    ϕ_exact, x_fine = get_exact_solution(u, x=x);
    fig = Figure()
    ax = Axis(fig[1,1], title="error", xlabel="x", ylabel="error")
    scatterlines!(ax, x, ϕ-ϕ_exact)
    ax2 = Axis(fig[1,2], title="ϕ", xlabel="x", ylabel="ϕ")
    scatterlines!(ax2, x, ϕ, label="TVD")
    lines!(ax2, x_fine, ϕ_exact, color=:black, linestyle=:dash, label="Exact")
    display(fig)
    @show x
    @show ϕ
end

##
# Create plots for each case
cases = ["i)", "i.2)", "ii)", "iii)"]
us = [0.1, 0.1,  2.5, 2.5]
Ns = [5, 8, 5, 90]

fig = Figure(size=(800, 600))
axes = [Axis(fig[i÷3+1, mod1(i, 2)], title = "Case $(cases[i]): u = $(us[i]), N = $(Ns[i])") for i in eachindex(cases)]

idxs = 2:8
limiters = [UD(),
            VanLeer(),
            VanAlbada(),
            Minmod(),
            Superbee(),
            Sweby(),
            UMIST(),
            QUICKlimiter()][idxs]
limiter_names = ["UD",
                 "VanLeer",
                 "VanAlbada",
                 "Minmod",
                 "Superbee",
                 "Sweby",
                 "UMIST",
                 "QUICK"][idxs]

for (i, u) in enumerate(us)
    @debug "case $i"
    ϕ_exact, x_fine = get_exact_solution(u)
    lines!(axes[i], x_fine, ϕ_exact, color=:black, linestyle=:dash, label="Exact")
    
    for (j, (limiter, name)) in enumerate(zip(limiters, limiter_names))
        discretization = TVD(limiter)
        @debug "Solving for $name"
        ϕ, x, its = solve_tvd(Ns[i], u, discretization)
        scatterlines!(axes[i], x, ϕ, label=name)
    end
end

# Add a single legend for all plots
fig[end+1,:] = Legend(fig, axes[1], "Methods", 
                framevisible=false, merge=true, 
                orientation=:horizontal, nbanks=2)
Label(fig[0,:], "Example 5.1")
# Adjust layout
for ax in axes
    ax.xlabel = "x"
    ax.ylabel = "ϕ"
end

fig
## save figure
save("plots/example_5-1_comparison.png", fig)
## Error analysis
us = [0.1, 2.5]
Ns = 2 .^ (3:8)
fig = Figure(size=(800, 600))
axes = [Axis(fig[i, 1],
            title = "u = $(us[i])",
            xscale=log2,
            yscale=log10,
            xlabel="N",
            ylabel="L2 error") for i in eachindex(us)]

conv_lines=[]
order_lines=[]
for (i, u) in enumerate(us)
    max_error = -Inf  # Track maximum error across all limiters
    for (j, limiter) in enumerate(limiters)
        L2_errors = []
        discretization = TVD(limiter)
        for N in Ns
            ϕ, x, its = solve_tvd(N, u, discretization);
            ϕ_exact, x_fine = get_exact_solution(u, x=x);
            dx = x[2] - x[1]
            if any(isnan.(ϕ))
                @error "ϕ is NaN for N = $N, limiter = $limiter"
            end
            push!(L2_errors, error(ϕ, ϕ_exact, dx))
        end
        max_error = max(max_error, maximum(L2_errors))  # Update maximum error
        l = lines!(axes[i], Ns, L2_errors, label=limiter_names[j])
        i == 1 && push!(conv_lines, l)
    end
    i == 2 && @show max_error
    # Add reference lines for convergence rates, scaled to current plot's data
    ref_x = [minimum(Ns), maximum(Ns)]
    ref_y_first = [max_error, max_error * (minimum(Ns)/maximum(Ns))]  # First-order reference line
    ref_y_second = [max_error, max_error * (minimum(Ns)/maximum(Ns))^2]  # Second-order reference line
    
    o1 = lines!(axes[i], ref_x, ref_y_first, linestyle=:dash, color=:black, label="First order")
    o2 = lines!(axes[i], ref_x, ref_y_second, linestyle=:dot, color=:black, label="Second order")
    i == 1 && push!(order_lines, o1)
    i == 1 && push!(order_lines, o2)
end
# Add a single legend for all plots
fig[end+1,:] = Legend(fig,
                    [conv_lines, order_lines], 
                    [limiter_names, ["First order", "Second order"]],
                    ["Limiter Functions", "Convergence"],
                    framevisible=false,
                    merge=true,
                    orientation=:vertical,
                    nbanks=8)
display(fig)
##
