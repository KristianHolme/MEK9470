function construct_linear_system(problem::LDCProblem, uvp::AbstractVector{T}) where T
    mesh = problem.mesh
    N = mesh.N
    ν = problem.ν
    ψ = problem.ψ
    bc_u = problem.bc_u
    bc_v = problem.bc_v
    ops = problem.ops

    # Handle AD compatibility: create new convection operator if types don't match
    conv_op = ConvectionOperator(mesh, T)
    update_convection_operators!(conv_op, mesh, uvp, ψ, bc_u, bc_v)
    Muu = conv_op.Muu
    Muv = conv_op.Muv
    Mvu = conv_op.Mvu
    Mvv = conv_op.Mvv
    
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

    #check if any matrix contains nans
    any(isnan.(Muu)) && @warn "Muu contains nans"
    any(isnan.(Muv)) && @warn "Muv contains nans"
    any(isnan.(Mvu)) && @warn "Mvu contains nans"
    any(isnan.(Mvv)) && @warn "Mvv contains nans"
    any(isnan.(lap_u)) && @warn "lap_u contains nans"
    any(isnan.(lap_v)) && @warn "lap_v contains nans"
    any(isnan.(cont_u)) && @warn "cont_u contains nans"
    any(isnan.(cont_v)) && @warn "cont_v contains nans"
    any(isnan.(dxp)) && @warn "dxp contains nans"
    any(isnan.(dyp)) && @warn "dyp contains nans"
    any(isnan.(lap_u_b)) && @warn "lap_u_b contains nans"
    any(isnan.(lap_v_b)) && @warn "lap_v_b contains nans"
    any(isnan.(cont_u_b)) && @warn "cont_u_b contains nans"
    any(isnan.(cont_v_b)) && @warn "cont_v_b contains nans"


    M = [Muu+ν*lap_u Muv dxp
        Mvu Mvv+ν*lap_v dyp
        cont_u cont_v spzeros(T, N^2, N^2)]
    # b = [lap_u_b + dxp_b
    # lap_v_b + dyp_b
    # cont_u_b + cont_v_b]
    b = Vector{T}([lap_u_b + dxp_b; lap_v_b + dyp_b; cont_u_b + cont_v_b])  # Convert to dense vector

    #replacing first cell continuity with pressure bc
    first_cell_continuity_index = 2 * N^2 + 1
    M[first_cell_continuity_index, :] .= 0
    b[first_cell_continuity_index] = 0
    dropzeros!(M)
    # dropzeros!(b)
    M[first_cell_continuity_index, first_cell_continuity_index] = 1

    return M, b
end