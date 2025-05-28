function construct_linear_system(problem::LDCProblem, uvp::AbstractVector)
    mesh = problem.mesh
    N = mesh.N
    ν = problem.ν
    ψ = problem.ψ
    bc_u = problem.bc_u
    bc_v = problem.bc_v
    ops = problem.ops

    update_convection_operators!(ops.convection, mesh, uvp, ψ, bc_u, bc_v)

    Muu = ops.convection.Muu
    Muv = ops.convection.Muv
    Mvu = ops.convection.Mvu
    Mvv = ops.convection.Mvv

    lap_u = ops.laplacian_u.M
    lap_u_b = ops.laplacian_u.b
    lap_v = ops.laplacian_v.M
    lap_v_b = ops.laplacian_v.b

    cont_u = ops.continuity_u.M
    cont_u_b = ops.continuity_u.b
    cont_v = ops.continuity_v.M
    cont_v_b = ops.continuity_v.b

    dxp = ops.dx_pressure.M
    dxp_b = ops.dx_pressure.b
    dyp = ops.dy_pressure.M
    dyp_b = ops.dy_pressure.b

    M = [Muu+ν*lap_u Muv dxp
        Mvu Mvv+ν*lap_v dyp
        cont_u cont_v spzeros(N^2, N^2)]
    b = [lap_u_b + dxp_b
        lap_v_b + dyp_b
        cont_u_b + cont_v_b]
    return M, b
end