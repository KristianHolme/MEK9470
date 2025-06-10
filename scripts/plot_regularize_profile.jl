using Dr
@quickactivate :MEK9470
using CairoMakie
##
set_theme!(theme_latexfonts())
##
x = 0:0.01:1 |> collect
y = similar(x)
@. y = 16 * (1 - x) .^ 2 .* x .^ 2

fig = Figure(size=(500, 400))
ax = Axis(fig[1, 1], xlabel="x", ylabel="u(x)", title="Regularized Lid Driven Cavity")

lines!(ax, x, y)
text!(ax, 0.26, 0.54, text="u(x) = 16(1 - x)² x²", fontsize=20)

fig
##
save("report2/plots/regularize_profile.svg", fig)