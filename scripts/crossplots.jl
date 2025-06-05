using DrWatson
@quickactivate :MEK9470
using CairoMakie
using Dates
using JLD2

"""
    create_crossplot(N::Int, Re::Int, limiters::Vector{<:AbstractLimiter}; 
                    quantity=:velocity_magnitude, figsize=(1200, 1200))

Create upper triangular crossplot comparing different limiters.

# Arguments
- `N`: Grid resolution
- `Re`: Reynolds number  
- `limiters`: Vector of limiters to compare
- `quantity`: What to plot (:velocity_magnitude, :u, :v, :vorticity, :streamfunction, :pressure)
- `figsize`: Figure size tuple

# Returns
- Makie Figure object

The plot shows:
- Diagonal: actual data using viridis colormap
- Off-diagonal: differences using vik colormap (white = zero)
"""
function create_crossplot(N::Int, Re::Int, limiters::Vector{<:AbstractLimiter};
    quantity=:velocity_magnitude, figsize=(1200, 1200))

    n_limiters = length(limiters)
    limiter_names = [string(lim) for lim in limiters]

    # Load solutions for all limiters
    solutions = Dict()
    meshes = Dict()

    println("Loading solutions for N=$N, Re=$Re...")
    for limiter in limiters
        dir = datadir("solutions", "N=$N", "Re=$Re")
        file_path = joinpath(dir, "$(string(limiter)).jld2")

        if !isfile(file_path)
            error("Solution file not found: $file_path")
        end

        data = load(file_path)
        solutions[limiter] = data["uvp"]
        mesh = CartesianMesh(N)
        meshes[limiter] = mesh
    end

    # Extract and compute the requested quantity for each limiter
    quantities = Dict()
    for limiter in limiters
        uvp = solutions[limiter]
        mesh = meshes[limiter]
        quantities[limiter] = extract_quantity(uvp, mesh, quantity)
    end

    # Create coordinate grids for plotting
    x = LinRange(0 + 1 / (2 * N), 1 - 1 / (2 * N), N)
    y = LinRange(0 + 1 / (2 * N), 1 - 1 / (2 * N), N)

    # Create figure with upper triangular layout
    fig = Figure(size=figsize)

    # Determine global colormap ranges
    all_values = vcat([vec(quantities[lim]) for lim in limiters]...)
    data_range = (minimum(all_values), maximum(all_values))

    # For differences, we want symmetric range around zero
    max_diff = 0.0
    for i in 1:n_limiters, j in i+1:n_limiters
        diff = quantities[limiters[i]] - quantities[limiters[j]]
        max_diff = max(max_diff, maximum(abs.(diff)))
    end
    diff_range = (-max_diff, max_diff)

    axes = Matrix{Any}(undef, n_limiters, n_limiters)

    # Create plots for upper triangle
    for i in 1:n_limiters
        for j in i:n_limiters
            # Calculate subplot position
            ax = Axis(fig[i, j],
                aspect=DataAspect(),
                title=i == j ? limiter_names[i] : "$(limiter_names[i]) - $(limiter_names[j])")

            axes[i, j] = ax

            if i == j
                # Diagonal: show actual data with viridis
                hm = heatmap!(ax, x, y, quantities[limiters[i]],
                    colormap=:viridis, colorrange=data_range)

                # Add colorbar for diagonal plots
                if j == n_limiters
                    Colorbar(fig[i, n_limiters+1], hm,
                        label=get_quantity_label(quantity))
                end
            else
                # Off-diagonal: show differences with vik (white = zero)
                diff = quantities[limiters[i]] - quantities[limiters[j]]
                hm = heatmap!(ax, x, y, diff,
                    colormap=:vik, colorrange=diff_range)

                # Add colorbar for difference plots (only rightmost column)
                if j == n_limiters && i == 1
                    Colorbar(fig[i, n_limiters+2], hm,
                        label="Difference")
                end
            end

            # Remove axes labels for cleaner look (except bottom row and left column)
            if i < n_limiters
                ax.xticklabelsvisible = false
            else
                ax.xlabel = "x"
            end

            if j > 1
                ax.yticklabelsvisible = false
            else
                ax.ylabel = "y"
            end
        end
    end

    # Add overall title
    quantity_str = get_quantity_title(quantity)
    Label(fig[0, :], "Limiter Comparison: $quantity_str (N=$N, Re=$Re)",
        fontsize=16, font="bold")

    # Add subtle grid lines to separate subplots
    for i in 1:n_limiters
        for j in i:n_limiters
            if haskey(axes, (i, j))
                lines!(axes[i, j], [0, 1, 1, 0, 0], [0, 0, 1, 1, 0],
                    color=:black, linewidth=1, alpha=0.3)
            end
        end
    end

    return fig
end

"""
    extract_quantity(uvp, mesh, quantity_symbol)

Extract the requested quantity from the solution vector.
"""
function extract_quantity(uvp, mesh, quantity_symbol)
    u, v, p = split_uvp(uvp, mesh)
    N = mesh.N

    U = reshape(u, N, N)
    V = reshape(v, N, N)
    P = reshape(p, N, N)

    if quantity_symbol == :velocity_magnitude
        return sqrt.(U .^ 2 + V .^ 2)
    elseif quantity_symbol == :u
        return U
    elseif quantity_symbol == :v
        return V
    elseif quantity_symbol == :pressure
        return P
    elseif quantity_symbol == :vorticity
        ω = compute_vorticity(uvp, mesh)
        return reshape(ω, N, N)
    elseif quantity_symbol == :streamfunction
        ψ = compute_streamfunction(uvp, mesh)
        return reshape(ψ, N, N)
    else
        error("Unknown quantity: $quantity_symbol")
    end
end

"""
Get appropriate label for quantity
"""
function get_quantity_label(quantity_symbol)
    labels = Dict(
        :velocity_magnitude => "|U|",
        :u => "u",
        :v => "v",
        :pressure => "p",
        :vorticity => "ω",
        :streamfunction => "ψ"
    )
    return get(labels, quantity_symbol, string(quantity_symbol))
end

"""
Get appropriate title for quantity
"""
function get_quantity_title(quantity_symbol)
    titles = Dict(
        :velocity_magnitude => "Velocity Magnitude",
        :u => "U-Velocity",
        :v => "V-Velocity",
        :pressure => "Pressure",
        :vorticity => "Vorticity",
        :streamfunction => "Streamfunction"
    )
    return get(titles, quantity_symbol, string(quantity_symbol))
end

"""
    batch_crossplots(N::Int, Re::Int, limiters::Vector; quantities=[:velocity_magnitude, :vorticity])

Generate crossplots for multiple quantities and save them.
"""
function batch_crossplots(N::Int, Re::Int, limiters::Vector;
    quantities=[:velocity_magnitude, :vorticity, :streamfunction])

    for quantity in quantities
        println("Creating crossplot for $quantity...")

        fig = create_crossplot(N, Re, limiters; quantity=quantity)

        # Save the plot
        filename = "crossplot_$(quantity)_N$(N)_Re$(Re).png"
        safesave(plotsdir("crossplots", filename), fig)

        println("  Saved: $(plotsdir("crossplots", filename))")
    end
end

# Example usage and test
if abspath(PROGRAM_FILE) == @__FILE__
    # Example parameters
    N = 128
    Re = 1000
    limiters = [VanLeer(), VanAlbada(), Minmod(), Superbee(), UMIST()]

    # Create individual crossplot
    fig = create_crossplot(N, Re, limiters; quantity=:velocity_magnitude)
    display(fig)

    # Create batch of crossplots
    # batch_crossplots(N, Re, limiters)
end

