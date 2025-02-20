using MEK9470
using CairoMakie

## Create range of r values to plot
r = collect(-1.0:0.1:8)

# Create figure
fig = Figure(size=(800, 600))
ax = Axis(fig[1,1], 
    title = "Flux Limiter Functions",
    xlabel = "r",
    ylabel = "ψ(r)")

# Plot each limiter function
limiters = [VanLeer(), VanAlbada(), Minmod(), Superbee(), Sweby(), UMIST(), QUICKlimiter()]
limiter_names = ["VanLeer", "VanAlbada", "Minmod", "Superbee", "Sweby", "UMIST", "QUICK"]

for (limiter, name) in zip(limiters, limiter_names)
    psi = [apply_limiter(limiter, ri) for ri in r]
    lines!(ax, r, psi, label=name)
end

# Add legend
axislegend(ax,position=:lt)

# Save figure
save("plots/limiter_functions.png", fig)

# Display figure
fig
