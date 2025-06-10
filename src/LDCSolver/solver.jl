function update_convection_operators!(op::ConvectionOperator,
    mesh::CartesianMesh,
    ϕu_f,
    ϕv_f,
    uvp,
    limiter::AbstractLimiter,
    bc_u::DomainBoundaryConditions,
    bc_v::DomainBoundaryConditions,
)
    # any(isnan.(uvp)) && @warn "uvp contains nans"
    N = mesh.N
    @assert length(uvp) == 3(N^2)
    u = @view uvp[1:N^2]
    v = @view uvp[N^2+1:2N^2]
    Muu, Muv = get_convection_operators(mesh, ϕu_f)
    Mvu, Mvv = get_convection_operators(mesh, ϕv_f)
    op.Muu = Muu
    op.Muv = Muv
    op.Mvu = Mvu
    op.Mvv = Mvv
    nothing
end

function construct_linear_system(problem::LDCProblem, uvp::AbstractVector{T}) where T
    mesh = problem.mesh
    N = mesh.N
    ν = problem.ν
    ψ = problem.ψ
    bc_u = problem.bc_u
    bc_v = problem.bc_v
    ops = problem.ops
    u, v, p = split_uvp(uvp, mesh)



    ϕu_f = compute_face_values(u, mesh, u, v, ψ, bc_u)
    ϕv_f = compute_face_values(v, mesh, u, v, ψ, bc_v)
    Muu, Muv = get_convection_operators(mesh, ϕu_f)
    Mvu, Mvv = get_convection_operators(mesh, ϕv_f)

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

    diagblock_u = Muu - lap_u .* ν
    diagblock_v = Mvv - lap_v .* ν

    mwi_face_corrections = get_face_correction_values(mesh, uvp, diagblock_u, diagblock_v, dxp, dyp)
    b_u_conv_corrections, b_v_conv_corrections = get_convection_correction_values(mwi_face_corrections, ϕu_f, ϕv_f, mesh)
    b_cont_correction = get_continuity_correction_values(mwi_face_corrections, mesh)


    #check if any matrix contains nans
    # any(isnan.(Muu)) && @warn "Muu contains nans"
    # any(isnan.(Muv)) && @warn "Muv contains nans"
    # any(isnan.(Mvu)) && @warn "Mvu contains nans"
    # any(isnan.(Mvv)) && @warn "Mvv contains nans"
    # any(isnan.(lap_u)) && @warn "lap_u contains nans"
    # any(isnan.(lap_v)) && @warn "lap_v contains nans"
    # any(isnan.(cont_u)) && @warn "cont_u contains nans"
    # any(isnan.(cont_v)) && @warn "cont_v contains nans"
    # any(isnan.(dxp)) && @warn "dxp contains nans"
    # any(isnan.(dyp)) && @warn "dyp contains nans"
    # any(isnan.(lap_u_b)) && @warn "lap_u_b contains nans"
    # any(isnan.(lap_v_b)) && @warn "lap_v_b contains nans"
    # any(isnan.(cont_u_b)) && @warn "cont_u_b contains nans"
    # any(isnan.(cont_v_b)) && @warn "cont_v_b contains nans"


    M = [diagblock_u Muv dxp
        Mvu diagblock_v dyp
        cont_u cont_v spzeros(T, N^2, N^2)]
    # b = [lap_u_b + dxp_b
    # lap_v_b + dyp_b
    # cont_u_b + cont_v_b]
    b = Vector{T}([
        ν * lap_u_b + dxp_b + b_u_conv_corrections;
        ν * lap_v_b + dyp_b + b_v_conv_corrections;
        cont_u_b + cont_v_b + b_cont_correction
    ])


    #replacing first cell continuity with pressure bc
    first_cell_continuity_index = 2 * N^2 + 1
    M[first_cell_continuity_index, :] .= T(0)
    b[first_cell_continuity_index] = T(0)
    dropzeros!(M)
    # dropzeros!(b)
    #average pressure zero
    # M[first_cell_continuity_index, first_cell_continuity_index:end] .= T(1)
    # first cell pressure zero
    M[first_cell_continuity_index, first_cell_continuity_index] = T(1)


    return M, b
end

function residuals(uvp, ldcprob)
    M, b = construct_linear_system(ldcprob, uvp)
    return M * uvp - b
end