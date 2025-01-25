using LinearAlgebra
using CairoMakie
N = 50 #number of control volumes
T_b = 100
T_inf = 20
L = 1.0 #m
n2 = 25 #h*P/(k*A)

function create_1D_grid(A, B, N)
    dx = (B - A)/N
    x = collect(LinRange(A +  dx/2, B - dx/2, N))
    return x, dx
end
x_A, x_B = 0.0, 0.0 + L
x, dx = create_1D_grid(x_A, x_B, N)

##
a_w = 1.0/dx .* ones(N)
a_e = 1.0/dx .* ones(N)
s_p = -n2*dx .* ones(N)
a_p = (a_w + a_e) .- s_p
s_u = n2*dx*T_inf.*ones(N)

#At boundary A
s_p[1] = -n2*dx - 2/dx
s_u[1] = n2*dx*T_inf + 2/dx * T_b
a_w[1] = 0
a_p[1] = a_w[1] + a_e[1] - s_p[1]

#At boundary B
s_p[end] = -n2*dx
s_u[end] = n2*dx*T_inf
a_e[end] = 0
a_p[end] = a_w[end] + a_e[end] - s_p[end]

M = Tridiagonal(-a_w[2:end], a_p, -a_e[1:end-1])
b = s_u

T = M\b

##
T_full = [T_b; T]
x_full = [x_A;x]

x_fine = collect(LinRange(x_A, x_B, 200))
T_exact = (T_b - T_inf) .* cosh.(sqrt(n2)*(L .- x_fine))/cosh(sqrt(n2)*L) .+ T_inf

fig = Figure()
ax = Axis(fig[1,1])
scatterlines!(ax, x_full, T_full)
lines!(ax, x_fine, T_exact, color=:red)
fig