using MEK9470
using CairoMakie

## Create range of r values to plot
r = collect(-1.0:0.1:8)
kwargs = (;xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5),
          xminorgridvisible=true, yminorgridvisible=true,
          xminorticksvisible=true, yminorticksvisible=true)
# Create figure
fig = Figure()
ax = Axis(fig[1,1], 
    title = "Flux Limiter Functions",
    xlabel = "r",
    ylabel = "ψ(r)"; kwargs...)

# Plot each limiter function
limiters = [VanLeer(), VanAlbada(), Minmod(), Superbee(), Sweby(1.5), UMIST(), QUICKlimiter()]
limiter_names = ["VanLeer", "VanAlbada", "Minmod", "Superbee", "Sweby(1.5)", "UMIST", "QUICK"]

for (limiter, name) in zip(limiters, limiter_names)
    psi = [apply_limiter(limiter, ri) for ri in r]
    lines!(ax, r, psi, label=name)
end

band!(ax, [0, 0.5, 1, 2, 8], [0, 0.5, 1, 1, 1], [0, 1, 1, 2, 2], color=:black, alpha=0.4,
    label="Second order region")

# Add legend
axislegend(ax,position=:rb, orientation=:horizontal, nbanks=3)
display(fig)
# Save figure
save("plots/limiter_functions.png", fig)
save("report/figures/limiter_functions.svg", fig)