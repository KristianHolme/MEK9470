function solve_up(ψ::AbstractLimiter, N::Int; Res=[1, 100, 250, 500, 1000, 2000, 5000])
    mesh = CartesianMesh(N)
    uvp = zeros(mesh.N^2 * 3) .+ 0.5
    termination_condition = RelNormSafeBestTerminationMode(norm)
    tol = 1e-4

    # Initialize continuation algorithm
    target_Re_list = sort(Int.(Res))  # Ensure sorted order and integer values
    last_successful_Re = 0
    current_Re = 1  # Always start from Re = 1
    completed_targets = Int[]
    failed_Re_list = Int[]  # Track failed Reynolds numbers
    failed_completely = false

    # Minimum Re gap to prevent infinite bisection (must be at least 1 for integers)
    min_Re_gap = 1

    # Adaptive stepping parameters
    max_step_factor = 0.5  # Maximum fraction to step toward next target
    safety_margin = 0.1    # Safety margin near failed Re values

    println("Starting continuation method for $(string(ψ)), N=$(N)")
    println("Target Reynolds numbers: $target_Re_list")

    while !failed_completely && length(completed_targets) < length(target_Re_list)
        ν = 1 / current_Re
        ldcprob = LDCProblem(mesh, ν, ψ)
        jac_sparsity = ADTypes.jacobian_sparsity(u -> residuals(u, ldcprob), uvp, TracerLocalSparsityDetector())
        ff = NonlinearFunction(residuals; jac_prototype=jac_sparsity)
        nlprob_sparse = NonlinearProblem(ff, uvp, ldcprob; abstol=tol, reltol=tol, termination_condition=termination_condition)

        println("  Attempting Re = $(current_Re)")
        t_start = time()
        sol = solve(nlprob_sparse, NewtonRaphson())
        t_end = time()
        t_elapsed = t_end - t_start

        if sol.retcode == ReturnCode.Success
            println("    ✓ Success at Re = $(current_Re)")
            uvp = sol.u  # Update for continuation
            last_successful_Re = current_Re

            # Check if we hit a target Reynolds number (exact match for integers)
            target_hit = false
            for target_Re in target_Re_list
                if target_Re ∉ completed_targets && current_Re == target_Re
                    # We hit a target - store solution and create plots
                    target_hit = true
                    completed_targets = sort([completed_targets; target_Re])

                    dir = datadir("solutions", "N=$(mesh.N)", "Re=$(target_Re)")
                    isdir(dir) || mkpath(dir)
                    res_norm = norm(residuals(uvp, ldcprob))
                    jldsave(joinpath(dir, "$(string(ψ)).jld2"); uvp, res_norm, t_elapsed)

                    # Create and save plot
                    fig = lid_driven_cavity_plot(sol.u, mesh, title="N=$(mesh.N), Re=$(target_Re), $(string(ψ))",
                        xlims=(0, 1), ylims=(0, 1), velocity_contours=true)
                    safesave(plotsdir("ldc", "N=$(mesh.N)", "Re=$(target_Re)", "$(string(ψ)).png"), fig)

                    println("    ★ Target Re = $(target_Re) completed and saved")
                    break
                end
            end

            # Determine next Re to attempt using adaptive stepping
            next_targets = filter(x -> x > current_Re && x ∉ completed_targets, target_Re_list)
            if !isempty(next_targets)
                next_target = minimum(next_targets)
                current_Re = compute_next_re(current_Re, next_target, failed_Re_list,
                    max_step_factor, safety_margin)
                println("    → Next Re = $(current_Re) (adaptive step toward $(next_target))")
            else
                # All targets completed
                break
            end

        else
            println("    ✗ Failed at Re = $(current_Re)")

            # Track failed Reynolds number
            if current_Re ∉ failed_Re_list
                failed_Re_list = sort([failed_Re_list; current_Re])
                println("    → Added Re = $(current_Re) to failed list: $(failed_Re_list)")
            end

            # Check if we can bisect
            if current_Re - last_successful_Re <= min_Re_gap
                println("    ! Re gap too small ($(current_Re - last_successful_Re) ≤ $(min_Re_gap)), abandoning further attempts")
                failed_completely = true
            else
                # Bisect between last successful and current Re (round to nearest integer)
                current_Re = round(Int, (last_successful_Re + current_Re) / 2)
                println("    → Bisecting, trying Re = $(current_Re)")
            end
        end
    end

    # Summary
    if length(completed_targets) == length(target_Re_list)
        println("✓ All target Reynolds numbers completed successfully!")
    else
        incomplete = filter(x -> x ∉ completed_targets, target_Re_list)
        println("⚠ Incomplete targets: $incomplete")
    end

    println("Completed targets: $(sort(completed_targets))")
    return completed_targets
end

"""
    compute_next_re(current_Re, target_Re, failed_Re_list, max_step_factor, safety_margin)

Compute the next Reynolds number to attempt using adaptive stepping.
This prevents jumping too aggressively toward targets after narrow successes.

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

    # Find the closest failed Re that's between current and target
    relevant_failures = filter(re -> current_Re < re < target_Re, failed_Re_list)

    if isempty(relevant_failures)
        # No failures between current and target - take a conservative step
        step_size = max(1, round(Int, gap_to_target * max_step_factor))
        next_Re = current_Re + step_size
    else
        # There are failures between current and target - be more careful
        closest_failure = minimum(relevant_failures)
        gap_to_failure = closest_failure - current_Re

        # Take a fraction of the gap to the closest failure, with safety margin
        safe_gap = gap_to_failure * (1 - safety_margin)
        step_size = max(1, round(Int, safe_gap))
        next_Re = current_Re + step_size

        # But don't exceed our normal stepping either
        max_normal_step = max(1, round(Int, gap_to_target * max_step_factor))
        next_Re = min(next_Re, current_Re + max_normal_step)
    end

    # Ensure we don't overshoot the target
    next_Re = min(next_Re, target_Re)

    # Ensure we make progress (at least +1)
    next_Re = max(next_Re, current_Re + 1)

    return next_Re
end