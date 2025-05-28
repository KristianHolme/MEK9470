using DrWatson
@quickactivate :MEK9470
##
mesh = CartesianMesh(4)

conv = ConvectionOperator(mesh)
uvp = 1:16*3 |> collect .|> Float64
update_convection_operators!(conv, mesh, uvp, VanLeer(), lid_driven_cavity_u_bc(1.0), lid_driven_cavity_v_bc())

ψ = VanLeer()
prob = LDCProblem(mesh, 1.0, ψ, 1.0)

M, b = construct_linear_system(prob, uvp)

x = M \ b

using NonlinearSolve

function residuals(uvp, prob)
    update_convection_operators!(prob.ops.convection, mesh, uvp, ψ, prob.bc_u, prob.bc_v)
    M, b = construct_linear_system(prob, uvp)
    return M * uvp - b
end

# function jacobian(uvp, prob)
#     u, v, p = split_uvp(uvp)
#     update_convection_operators!(prob.ops.convection, mesh, uvp, ψ, prob.bc_u, prob.bc_v)
#     M, b = construct_linear_system(prob, uvp)
#     return M
# end

uvp = rand(16 * 3)
nlprob = NonlinearProblem(residuals, uvp, prob)
sol = solve(nlprob)