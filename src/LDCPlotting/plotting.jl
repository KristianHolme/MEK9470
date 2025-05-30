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
    title="Lid Driven Cavity", show_streamlines=true,
    figsize=(800, 400))
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
        aspect=DataAspect())

    # Background: velocity magnitude
    hm1 = heatmap!(ax1, x, y, vel_mag, colormap=:viridis)
    Colorbar(fig[1, 2], hm1, label="Velocity Magnitude")

    # Overlay: velocity vectors (subsample for clarity)
    skip = max(1, N ÷ 8)  # Show every nth vector
    arrows!(ax1, x[1:skip:end], y[1:skip:end],
        U[1:skip:end, 1:skip:end], V[1:skip:end, 1:skip:end],
        arrowsize=10, lengthscale=0.3, color=:white, linewidth=1.5)

    # Optional streamlines
    if show_streamlines
        contour!(ax1, x, y, vel_mag, levels=10, color=:white, alpha=0.6, linewidth=0.8)
    end

    # Pressure plot
    ax2 = Axis(fig[1, 3],
        title="Pressure Field",
        xlabel="x", ylabel="y",
        aspect=DataAspect())

    hm2 = heatmap!(ax2, x, y, P, colormap=:viridis)
    Colorbar(fig[1, 4], hm2, label="Pressure")

    # Add cavity boundary visualization
    for ax in [ax1, ax2]
        # Draw cavity walls
        lines!(ax, [0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color=:black, linewidth=3)
        # Highlight moving lid
        lines!(ax, [0, 1], [1, 1], color=:red, linewidth=4)
        xlims!(ax, 0, 1)
        ylims!(ax, 0, 1)
    end

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
    plot_streamfunction(uvp, mesh)

Compute and plot the streamfunction for the lid driven cavity.
The streamfunction ψ satisfies: ∂ψ/∂y = u, ∂ψ/∂x = -v
"""
function plot_streamfunction(uvp::AbstractVector, mesh::CartesianMesh)
    u, v, p = split_uvp(uvp, mesh)
    N = mesh.N

    # Create coordinate grids at cell centers
    x = LinRange(0 + 1 / (2 * N), 1 - 1 / (2 * N), N)
    y = LinRange(0 + 1 / (2 * N), 1 - 1 / (2 * N), N)

    U = reshape(u, N, N)
    V = reshape(v, N, N)

    # Compute streamfunction by integrating velocity field
    # ψ = ∫u dy (integrating along y direction)
    dx = x[2] - x[1]
    dy = y[2] - y[1]

    ψ = zeros(N, N)
    # Integrate u in y-direction  
    for i in 1:N
        for j in 2:N
            ψ[i, j] = ψ[i, j-1] + U[i, j-1] * dy
        end
    end

    fig = Figure(size=(600, 500))
    ax = Axis(fig[1, 1],
        title="Streamfunction",
        xlabel="x", ylabel="y",
        aspect=DataAspect())

    # Plot streamlines
    contour!(ax, x, y, ψ', levels=20, color=:black)
    contourf!(ax, x, y, ψ', levels=20, colormap=:RdBu_r)

    # Add cavity boundary
    lines!(ax, [0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color=:black, linewidth=3)
    lines!(ax, [0, 1], [1, 1], color=:red, linewidth=4)  # Moving lid

    xlims!(ax, 0, 1)
    ylims!(ax, 0, 1)

    Colorbar(fig[1, 2], label="Streamfunction ψ")

    return fig
end