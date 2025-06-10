using MEK9470
using CairoMakie
set_theme!(theme_latexfonts())
## Create range of r values to plot
r = collect(-1.0:0.001:8)
kwargs = (; xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5),
    xminorgridvisible=true, yminorgridvisible=true,
    xminorticksvisible=true, yminorticksvisible=true)
# Create figure
fig = Figure(size=(500, 400))
ax = Axis(fig[1, 1], limits=((0, 3), nothing),
    title="Flux Limiter Functions",
    xlabel="r",
    ylabel="ψ(r)"; kwargs...)

# Plot each limiter function
limiters = [VanLeer(), VanAlbada(), Minmod(), Superbee(), Sweby(1.5), UMIST(), QUICKlimiter()]
limiter_names = ["VanLeer", "VanAlbada", "Minmod", "Superbee", "Sweby(1.5)", "UMIST", "QUICK"]

band!(ax, [0, 0.5, 1, 2, 8], [0, 0.5, 1, 1, 1], [0, 1, 1, 2, 2], color=:black, alpha=0.4,
    label="Second order region")
for (limiter, name) in zip(limiters, limiter_names)
    psi = [apply_limiter(limiter, ri) for ri in r]
    lines!(ax, r, psi, label=name)
end


# Add legend
axislegend(ax, position=:rb, orientation=:horizontal, nbanks=3)
display(fig)
## Save figure
safesave("plots/limiter_functions.svg", fig)
# safesave("report/figures/limiter_functions.svg", fig)
safesave("report2/plots/limiter_functions.svg", fig)
