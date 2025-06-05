function simple_uvp_plot(uvp::AbstractVector, mesh::CartesianMesh)
    u, v, p = split_uvp(uvp, mesh)
    N = mesh.N
    fig = Figure()
    ax_u = Axis(fig[1, 1], title="u")
    ax_v = Axis(fig[1, 2], title="v")
    ax_p = Axis(fig[2, 1], title="p")
    heatmap!(ax_u, reshape(u, N, N))
    heatmap!(ax_v, reshape(v, N, N))
    heatmap!(ax_p, reshape(p, N, N))
    fig
end

"""
    lid_driven_cavity_plot(uvp, mesh; title="Lid Driven Cavity", show_streamlines=true)

Create a comprehensive plot of the lid driven cavity solution showing:
- Velocity magnitude as background
- Velocity vectors as arrows  
- Optional streamlines
- Pressure field in separate subplot
"""
function lid_driven_cavity_plot(uvp::AbstractVector, mesh::CartesianMesh;
    title="Lid Driven Cavity", velocity_contours=true,
    figsize=(800, 800), pressure_lims=nothing, xlims=nothing, ylims=nothing)
    u, v, p = split_uvp(uvp, mesh)
    N = mesh.N

    # Create coordinate grids
    x = LinRange(0 + 1 / (2 * N), 1 - 1 / (2 * N), N)
    y = LinRange(0 + 1 / (2 * N), 1 - 1 / (2 * N), N)

    # Reshape velocity and pressure fields (same as simple_uvp_plot)
    U = reshape(u, N, N)
    V = reshape(v, N, N)
    P = reshape(p, N, N)

    # Compute velocity magnitude
    vel_mag = sqrt.(U .^ 2 + V .^ 2)

    fig = Figure(size=figsize)

    # Main velocity plot
    ax1 = Axis(fig[1, 1],
        title="Velocity Field",
        xlabel="x", ylabel="y",
        aspect=DataAspect(), limits=(xlims, ylims))

    # Background: velocity magnitude
    hm1 = heatmap!(ax1, x, y, vel_mag, colormap=:viridis)
    Colorbar(fig[1, 2], hm1, label="Velocity Magnitude")

    # Overlay: velocity vectors (subsample for clarity)
    skip = max(1, N ÷ 8)  # Show every nth vector
    arrows!(ax1, x[1:skip:end], y[1:skip:end],
        U[1:skip:end, 1:skip:end], V[1:skip:end, 1:skip:end],
        arrowsize=10, lengthscale=0.3, color=:white, linewidth=1.5)

    # Optional streamlines
    if velocity_contours
        contour!(ax1, x, y, vel_mag, levels=10, color=:white, alpha=0.6, linewidth=0.8)
    end

    # Pressure plot
    ax2 = Axis(fig[1, 3],
        title="Pressure Field",
        xlabel="x", ylabel="y",
        aspect=DataAspect(), limits=(xlims, ylims))

    if pressure_lims === nothing
        hm2 = heatmap!(ax2, x, y, P, colormap=:viridis)
    else
        hm2 = heatmap!(ax2, x, y, P, colormap=:viridis, colorrange=pressure_lims)
    end
    Colorbar(fig[1, 4], hm2, label="Pressure")

    # Streamfunction plot
    ψ = compute_streamfunction(uvp, mesh)
    Ψ = reshape(ψ, N, N)
    ax3 = Axis(fig[2, 1], title="Streamfunction ψ", xlabel="x", ylabel="y", aspect=DataAspect())
    # psi_plot = contourf!(ax3, x, y, Ψ, levels=20, colormap=:balance)
    psi_plot = heatmap!(ax3, x, y, Ψ, colormap=:viridis)
    contour!(ax3, x, y, Ψ, levels=[0.0, -10, -200, -500, -1000, -1300], color=:black, linewidth=1)
    # arrows!(ax3, x[1:skip:end], y[1:skip:end],
    #     U[1:skip:end, 1:skip:end], V[1:skip:end, 1:skip:end],
    #     arrowsize=10, lengthscale=0.3, color=:red, linewidth=1.5)
    Colorbar(fig[2, 2], psi_plot, label="ψ")

    # Vorticity plot
    ω = compute_vorticity(uvp, mesh)
    Ω = reshape(ω, N, N)
    ax4 = Axis(fig[2, 3], title="Vorticity ω", xlabel="x", ylabel="y", aspect=DataAspect())
    omega_plot = contourf!(ax4, x, y, Ω, levels=20, colormap=:viridis)
    # omega_plot = heatmap!(ax4, x, y, Ω, colormap=:viridis)

    contour!(ax4, x, y, Ω, levels=10, color=:black, linewidth=1)
    Colorbar(fig[2, 4], omega_plot, label="ω")

    # # Add cavity boundary visualization
    # for ax in [ax1, ax2]
    #     # Draw cavity walls
    #     lines!(ax, [0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color=:black, linewidth=3)
    #     # Highlight moving lid
    #     lines!(ax, [0, 1], [1, 1], color=:red, linewidth=4)
    #     xlims!(ax, 0, 1)
    #     ylims!(ax, 0, 1)
    # end

    Label(fig[0, :], title, fontsize=16, font="bold")

    return fig
end

"""
    plot_velocity_profiles(uvp, mesh; y_positions=[0.5], x_positions=[0.5])

Plot velocity profiles at specified horizontal and vertical lines for validation.
Useful for comparing with benchmark data.
"""
function plot_velocity_profiles(uvp::AbstractVector, mesh::CartesianMesh;
    y_positions=[0.5], x_positions=[0.5])
    u, v, p = split_uvp(uvp, mesh)
    N = mesh.N

    # Create coordinate grids at cell centers
    x = LinRange(0 + 1 / (2 * N), 1 - 1 / (2 * N), N)
    y = LinRange(0 + 1 / (2 * N), 1 - 1 / (2 * N), N)

    U = reshape(u, N, N)
    V = reshape(v, N, N)

    fig = Figure(size=(800, 400))

    # U-velocity profiles at horizontal lines
    ax1 = Axis(fig[1, 1],
        title="u-velocity profiles",
        xlabel="x", ylabel="u")

    for y_pos in y_positions
        j = argmin(abs.(y .- y_pos))
        u_profile = U[:, j]
        lines!(ax1, x, u_profile, label="y = $(round(y[j], digits=2))")
    end
    axislegend(ax1)

    # V-velocity profiles at vertical lines  
    ax2 = Axis(fig[1, 2],
        title="v-velocity profiles",
        xlabel="y", ylabel="v")

    for x_pos in x_positions
        i = argmin(abs.(x .- x_pos))
        v_profile = V[i, :]
        lines!(ax2, y, v_profile, label="x = $(round(x[i], digits=2))")
    end
    axislegend(ax2)

    return fig
end

"""
    plot_streamfunction(uvp, mesh; method=:simpson, levels=20)

Compute and plot the streamfunction for the lid driven cavity using advanced integration methods.
The streamfunction ψ satisfies: ∂ψ/∂x = v, ∂ψ/∂y = -u
"""
function plot_streamfunction(uvp::AbstractVector, mesh::CartesianMesh; method::Symbol=:simpson, levels::Int=20)
    N = mesh.N

    # Create coordinate grids at cell centers
    x = LinRange(0 + 1 / (2 * N), 1 - 1 / (2 * N), N)
    y = LinRange(0 + 1 / (2 * N), 1 - 1 / (2 * N), N)

    # Use the improved streamfunction computation
    ψ_vec = compute_streamfunction(uvp, mesh; method=method)
    ψ = reshape(ψ_vec, mesh)

    fig = Figure(size=(600, 500))
    ax = Axis(fig[1, 1],
        title="Streamfunction ($(string(method)) integration)",
        xlabel="x", ylabel="y",
        aspect=DataAspect())

    # Plot streamlines
    contour!(ax, x, y, ψ', levels=levels, color=:black, linewidth=1.5)
    contourf!(ax, x, y, ψ', levels=levels, colormap=:viridis)

    # Add cavity boundary
    lines!(ax, [0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color=:black, linewidth=3)
    lines!(ax, [0, 1], [1, 1], color=:red, linewidth=4)  # Moving lid

    xlims!(ax, 0, 1)
    ylims!(ax, 0, 1)

    Colorbar(fig[1, 2], label="Streamfunction ψ")

    return fig
end

"""
    plot_flow_analysis(uvp, mesh; method=:simpson)

Create a comprehensive flow analysis plot showing streamfunction, velocity potential, 
and vorticity fields.
"""
function plot_flow_analysis(uvp::AbstractVector, mesh::CartesianMesh;)
    N = mesh.N

    # Create coordinate grids
    x = LinRange(0 + 1 / (2 * N), 1 - 1 / (2 * N), N)
    y = LinRange(0 + 1 / (2 * N), 1 - 1 / (2 * N), N)

    # Compute flow functions
    ψ_vec = compute_streamfunction(uvp, mesh)
    ω_vec = compute_vorticity(uvp, mesh)

    # Reshape for plotting
    Ψ = reshape(ψ_vec, N, N)
    Ω = reshape(ω_vec, N, N)

    fig = Figure(size=(1200, 400))

    # Streamfunction plot
    ax1 = Axis(fig[1, 1], title="Streamfunction ψ", xlabel="x", ylabel="y", aspect=DataAspect())
    contourf!(ax1, x, y, Ψ, levels=20, colormap=:viridis)
    contour!(ax1, x, y, Ψ, levels=15, color=:black, linewidth=1)
    lines!(ax1, [0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color=:black, linewidth=2)
    Colorbar(fig[1, 2], label="ψ")

    # Vorticity plot
    ax3 = Axis(fig[1, 5], title="Vorticity ω", xlabel="x", ylabel="y", aspect=DataAspect())
    ω_max = maximum(abs.(Ω))
    contourf!(ax3, x, y, Ω, levels=20, colormap=:viridis)
    contour!(ax3, x, y, Ω, levels=10, color=:black, linewidth=1)
    lines!(ax3, [0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color=:black, linewidth=2)
    Colorbar(fig[1, 6], label="ω")

    for ax in [ax1, ax2, ax3]
        xlims!(ax, 0, 1)
        ylims!(ax, 0, 1)
    end

    return fig
end