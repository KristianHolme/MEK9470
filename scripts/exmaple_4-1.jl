using LinearAlgebra
using CairoMakie
N = 5 #number of control volumes
k = 1000 # W/mK
T_a = 100
T_b = 500
A = 10e-3 # m^2
L = 0.5 #m

function create_1D_grid(A, B, N)
    dx = (B - A)/N
    x = collect(LinRange(A +  dx/2, B - dx/2, N))
    return x, dx
end
A, B = 0.0, 0.5
x, dx = create_1D_grid(A, B, 5)

##
a_w = k * A / dx .* ones(N)
a_e = k * A / dx .* ones(N)
a_p = (a_w + a_e) .* ones(N)
s_p = zeros(N)
s_u = zeros(N)

#At boundary A
s_p[1] = -2*k*A/dx
s_u[1] = 2*k*A/dx * T_a
a_w[1] = 0
a_p[1] = a_e[1] - s_p[1]

#At boundary B
s_p[end] = -2*k*A/dx
s_u[end] = 2*k*A/dx*T_b
a_e[end] = 0
a_p[end] = a_w[end] + a_e[end] - s_p[end]

M = Tridiagonal(-a_w[2:end], a_p, -a_e[1:end-1])
b = s_u

T = M\b

##
T_full = [T_a; T; T_b]
x_full = [A;x;B]

scatterlines(x_full, T_full)