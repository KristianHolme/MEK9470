function solve_up(ψ::AbstractLimiter, N::Int; Res=[1, 10, 20, 40, 70, 100, 150, 200, 250, 350, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000],
    save=true)
    mesh = CartesianMesh(N)
    uvp = zeros(mesh.N^2 * 3) .+ 0.5
    tol = 1e-4

    failed = false
    println("Starting continuation method for $(string(ψ)), N=$(N)")
    println("Reynolds numbers: $Res")

    for Re in Res
        ν = 1 / Re
        ldcprob = LDCProblem(mesh, ν, ψ)
        jac_sparsity = ADTypes.jacobian_sparsity(u -> residuals(u, ldcprob), uvp, TracerLocalSparsityDetector())
        ff = NonlinearFunction(residuals; jac_prototype=jac_sparsity)
        termination_condition = RelNormSafeBestTerminationMode(norm)
        nlprob_sparse = NonlinearProblem(ff, uvp, ldcprob; abstol=tol, reltol=tol, termination_condition=termination_condition)

        println("  Attempting Re = $(Re)")
        t_start = time()
        sol = solve(nlprob_sparse, NewtonRaphson())
        t_end = time()
        t_elapsed = t_end - t_start

        if sol.retcode == ReturnCode.Success
            println("    ✓ Success at Re = $(Re)")
            uvp = sol.u  # Update for continuation
            if save
                dir = datadir("solutions", "N=$(mesh.N)", "Re=$(Re)")
                isdir(dir) || mkpath(dir)
                res_norm = norm(residuals(uvp, ldcprob))
                jldsave(joinpath(dir, "$(string(ψ)).jld2"); uvp, res_norm, t_elapsed)

                # Create and save plot
                fig = lid_driven_cavity_plot(sol.u, mesh, title="N=$(mesh.N), Re=$(Re), $(string(ψ))",
                    xlims=(0, 1), ylims=(0, 1), velocity_contours=true)
                safesave(plotsdir("ldc", "N=$(mesh.N)", "Re=$(Re)", "$(string(ψ)).png"), fig)

                stream = compute_streamfunction(sol.u, mesh)
                fig = Figure(size=(600, 600))
                ax = Axis(fig[1, 1], title="Streamfunction, Re = $(Re), ψ = $(string(ψ))", xlabel="x", ylabel="y", aspect=DataAspect())
                stream_mat = reshape(stream, mesh.N, mesh.N)
                coords = cell_center.(1:N, 1:N, Ref(mesh))
                xs = getindex.(coords, Ref(1))
                ys = getindex.(coords, Ref(2))
                hm = heatmap!(ax, xs, ys, stream_mat, colormap=:coolwarm)
                I = argmin(stream_mat)
                i, j = I[1], I[2]
                ccx, ccy = cell_center(i, j, mesh)
                sc = scatter!(ax, [ccx], [ccy], color=:magenta, markersize=10, label="Minimum")
                text!(ax, ccx + 0.01, ccy + 0.01, text=string("($(round(ccx, digits=3)), $(round(ccy, digits=3)))"), fontsize=12)
                Colorbar(fig[1, 2], hm)
                safesave(plotsdir("ldc", "N=$(mesh.N)", "Re=$(Re)", "$(string(ψ))_streamfunction.png"), fig)
            end
        else
            println("    ✗ Failed at Re = $(Re), retcode: $(sol.retcode)")
            failed = true
            break
        end
    end
end

"""
    compute_next_re(current_Re, target_Re, failed_Re_list, max_step_factor, safety_margin)

Compute the next Reynolds number to attempt using efficient stepping with minimum step size.
- Minimum step = 10, except if within 10 of target, then step directly to target
- This prevents wasteful small steps near targets

# Arguments
- `current_Re`: Current successful Reynolds number
- `target_Re`: Next target Reynolds number to reach
- `failed_Re_list`: List of previously failed Reynolds numbers
- `max_step_factor`: Maximum fraction of the gap to step (0.0-1.0)
- `safety_margin`: Safety margin as fraction when near failed Re values

# Returns
- Next Reynolds number to attempt (integer)
"""
function compute_next_re(current_Re::Int, target_Re::Int, failed_Re_list::Vector{Int},
    max_step_factor::Float64, safety_margin::Float64)

    # Calculate the gap to the target
    gap_to_target = target_Re - current_Re

    # Minimum step size (10) unless close to target
    min_step = 10

    # If within min_step of target, go directly to target
    if gap_to_target <= min_step
        return target_Re
    end

    # Find the closest failed Re that's between current and target
    relevant_failures = filter(re -> current_Re < re < target_Re, failed_Re_list)

    if isempty(relevant_failures)
        # No failures between current and target - take a step based on gap
        step_size = max(min_step, round(Int, gap_to_target * max_step_factor))
        next_Re = current_Re + step_size
    else
        # There are failures between current and target - be more careful
        closest_failure = minimum(relevant_failures)
        gap_to_failure = closest_failure - current_Re

        # If the closest failure is very close, jump to it
        if gap_to_failure <= min_step
            step_size = gap_to_failure
        else
            # Use normal minimum step size, but respect failure boundaries
            step_size = max(min_step, round(Int, gap_to_target * max_step_factor))
        end
        next_Re = current_Re + step_size
    end

    # Ensure we don't overshoot the target
    next_Re = min(next_Re, target_Re)

    # Ensure we make progress (at least +1)
    next_Re = max(next_Re, current_Re + 1)

    return next_Re
end