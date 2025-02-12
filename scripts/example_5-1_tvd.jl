using MEK9470
using LinearAlgebra
using CairoMakie
##Parameters
A = 0.0
B = 1.0
ρ = 1
Γ = 0.1
ϕ_A = 1.0
ϕ_B = 0.0
F = ρ*u
##
N = 5
x, dx = create_1D_grid(A, B, N)

D = Γ / dx

a_w = D + F
a_e = D



function get_exact_solution(u)
    F = ρ*u
    x_fine = collect(LinRange(A, B, 100))
    ϕ_exact = (ϕ_B - ϕ_A) *(exp.(F.*x_fine./Γ) .- 1) ./ (exp.(F*B./Γ) .- 1) .+ ϕ_A
    return x_fine, ϕ_exact
end