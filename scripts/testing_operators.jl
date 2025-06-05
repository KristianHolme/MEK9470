using DrWatson
@quickactivate :MEK9470
using CairoMakie
##
mesh = CartesianMesh(8)

dxp, b = get_pressure_dx_operator(mesh)
##
fig = Figure()
ax = Axis(fig[1, 1])
sp = spy!(ax, dxp[1:10, 1:10])
Colorbar(fig[1, 2], sp)
fig
##
dyp, b = get_pressure_dy_operator(mesh)
##
fig = Figure()
ax = Axis(fig[1, 1])
sp = spy!(ax, dyp)
Colorbar(fig[1, 2], sp)
fig
##
cont_x, b = get_continuity_operator_x(mesh)
##
fig = Figure()
ax = Axis(fig[1, 1])
sp = spy!(ax, cont_x)
Colorbar(fig[1, 2], sp)
fig
##
cont_y = get_continuity_operator_y(mesh)
##
fig = Figure()

lap, b = get_laplacian_operator(mesh, 0)
##
fig = Figure()
ax = Axis(fig[1, 1])
sp = spy!(ax, lap)
Colorbar(fig[1, 2], sp)
fig
##

