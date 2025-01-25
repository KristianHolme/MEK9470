using LinearAlgebra
using CairoMakie
N = 5 #number of control volumes
k = 0.5 # W/mK
q = 1000e3 # kW/m^3
T_a = 100
T_b = 200
A = 1 # m^2
L = 0.02 #m

function create_1D_grid(A, B, N)
    dx = (B - A)/N
    x = collect(LinRange(A +  dx/2, B - dx/2, N))
    return x, dx
end
x_A, x_B = 0.0, 0.02
x, dx = create_1D_grid(x_A, x_B, 5)

##
a_w = k * A / dx .* ones(N)
a_e = k * A / dx .* ones(N)
a_p = (a_w + a_e) .* ones(N)
s_p = zeros(N)
s_u = q*A*dx.*ones(N)

#At boundary A
s_p[1] = -2*k*A/dx
s_u[1] = 2*k*A/dx * T_a + q*A*dx
a_w[1] = 0
a_p[1] = a_w[1] + a_e[1] - s_p[1]

#At boundary B
s_p[end] = -2*k*A/dx
s_u[end] = 2*k*A/dx*T_b + q*A*dx
a_e[end] = 0
a_p[end] = a_w[end] + a_e[end] - s_p[end]

M = Tridiagonal(-a_w[2:end], a_p, -a_e[1:end-1])
b = s_u

T = M\b

##
T_full = [T_a; T; T_b]
x_full = [x_A;x;x_B]

x_fine = collect(LinRange(x_A, x_B, 200))
T_exact = ((T_b - T_a)/L .+ q/(2*k).*(L .- x_fine)).*x_fine .+ T_a

fig = Figure()
ax = Axis(fig[1,1])
scatterlines!(ax, x_full, T_full)
lines!(ax, x_fine, T_exact, color=:red)
fig