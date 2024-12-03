using LinearAlgebra, FFMPEG, BenchmarkTools, PyPlot, Optimization, Plots
using Base.MathConstants: π

# Constants and parameters for the shaping functions
const D = 750.0       # Characteristic downburst diameter
const zm = 0.025 * D  # Height of maximum horizontal velocity
const rm = 1.125 * D  # Radial distance of maximum horizontal velocity
const λ = 0.2  # Scaling factor (to be calibrated)

# Shaping function parameters
const ζ1 = 0.133
const ζ2 = 1.1534
const α = 1.8
const β = 300; 

function WindField(pos::Vector{Float64})
    E, N, z = pos[1:3]
    r = sqrt.(E .^ 2 + N .^ 2)  # radial distance from core
    # eq. A.11
    u = (λ*r / 2) *
        (
            exp(ζ1 * (z / zm)) - exp(ζ2 * (z / zm))
        ) *
        (
            exp(
            (
                1 - ((r^2) / (β^2))^α
            ) / α
        )
        )
    #eq. A.12
    w = (
            (zm / ζ1) * (exp(ζ1 * (z / zm)) - 1)
            -
            (zm / ζ2) * (exp(ζ2 * (z / zm) - 1))
        ) *
        (
            1 - (r^2 / β^2)^α
        ) *
        exp(
            (1 - (r^2 / β^2)^α)
            /
            α
        )

    # cylindrical to ac
    if (E != 0 || N != 0)
        θ = atan(N,E)
    else
        θ = 45*π/180
    end

    return [u*cos(θ), u*sin(θ), w]
end

# generic wind field 
function GenericWindField(pos)
    E, N, U = pos[1:3]
    W_N = U / 10
    W_E = U / 10
    W_U = 0
    return [W_E; W_N; W_U]
end
# EU quiver wind plot
function plot_wind_direction_VF(domain_x, domain_z, resolution=50, arrow_scale=0.2)
    # Generate x and z vectors
    x_vals = LinRange(domain_x[1], domain_x[2], resolution)
    z_vals = LinRange(domain_z[1], domain_z[2], resolution)
    
    # Create meshgrid for x and z
    x_coords = repeat(x_vals, 1, resolution)'
    z_coords = repeat(z_vals', resolution, 1)

    # Initialize wind field components
    u = zeros(resolution, resolution)
    w = zeros(resolution, resolution)

    # Compute wind field values at each grid point
    for i in 1:resolution
        for j in 1:resolution
            u[i, j], _, w[i, j] = WindField([x_coords[i, j], 0.0, z_coords[i, j]])
        end
    end

    # Calculate magnitude of the wind field
    mag = sqrt.(u .^ 2 + w .^ 2)

    # Normalize vectors for quiver plot and scale them
    u_scaled = arrow_scale .* (u ./ (mag .+ eps()))  # Avoid division by zero
    w_scaled = arrow_scale .* (w ./ (mag .+ eps()))

    # Set up plot
    gr(legend=false, dpi=600)

    # Heatmap for magnitude
    heatmap(
        x_vals, z_vals, mag',
        color=:viridis,
        title="Wind Field (Vertical Slice at N=0)",
        xlabel="East [m]",
        ylabel="Up [m]"
    )

    # Overlay quiver plot with scaled vectors
    quiver!(
        x_coords, z_coords, quiver=(u_scaled, w_scaled), color=:black
    )

    # Save and display the plot
    Plots.savefig("./figures/WindVF.png")
    show()
    return nothing
end
function plot_wind_top_down(domain_x, domain_y, height=50.0, resolution=50, arrow_scale=0.2)
    # Generate x and y vectors
    x_vals = LinRange(domain_x[1], domain_x[2], resolution)
    y_vals = LinRange(domain_y[1], domain_y[2], resolution)
    # Meshgrid
    x_coords = repeat(x_vals, 1, resolution)
    y_coords = repeat(y_vals', resolution, 1)

    # Compute wind field at specified height
    u = zeros(resolution, resolution)
    v = zeros(resolution, resolution)
    for i in 1:resolution
        for j in 1:resolution
            u[i, j], v[i, j], _ = WindField([x_coords[i, j], y_coords[i, j], height])
        end
    end

    # Calculate wind magnitude
    mag = sqrt.(u .^ 2 + v .^ 2)

    # Normalize the vectors and scale them for the quiver plot
    u_scaled = arrow_scale .* (u ./ (mag .+ eps()))  
    v_scaled = arrow_scale .* (v ./ (mag .+ eps()))

    # Set up plot
    gr(legend=false, dpi=600)

    # Heatmap for wind magnitude
    heatmap(
        x_vals, y_vals, mag,
        color=:viridis,
        title="Wind Field at Height = $height m",
        xlabel="East [m]",
        ylabel="North [m]"
    )

    # Overlay quiver plot with scaled vectors
    quiver!(x_coords, y_coords, quiver=(u_scaled, v_scaled), color=:black)

    # Save and display plot
    Plots.savefig("./figures/WindTopDown.png")
    show()
    return nothing
end

function plot_wind_top_down_with_path(domain_x, domain_y, history::Matrix{Float64}, height=50.0, resolution=50, arrow_scale=0.2)
    # Generate x and y vectors for the vector field
    x_vals = LinRange(domain_x[1], domain_x[2], resolution)
    y_vals = LinRange(domain_y[1], domain_y[2], resolution)
    # Meshgrid for vector field
    x_coords = repeat(x_vals, 1, resolution)
    y_coords = repeat(y_vals', resolution, 1)

    # Compute wind field at specified height
    u = zeros(resolution, resolution)
    v = zeros(resolution, resolution)
    for i in 1:resolution
        for j in 1:resolution
            u[i, j], v[i, j], _ = WindField([x_coords[i, j], y_coords[i, j], height])
        end
    end

    # Calculate wind magnitude
    mag = sqrt.(u .^ 2 + v .^ 2)

    # Normalize the vectors and scale them for the quiver plot
    u_scaled = arrow_scale .* (u ./ (mag .+ eps()))  # Avoid division by zero
    v_scaled = arrow_scale .* (v ./ (mag .+ eps()))

    # Generate high-resolution grid for heatmap
    heatmap_resolution = resolution * 10
    heatmap_x_vals = LinRange(domain_x[1], domain_x[2], heatmap_resolution)
    heatmap_y_vals = LinRange(domain_y[1], domain_y[2], heatmap_resolution)
    heatmap_mag = zeros(heatmap_resolution, heatmap_resolution)

    for i in 1:heatmap_resolution
        for j in 1:heatmap_resolution
            heatmap_mag[i, j], _, _ = WindField([heatmap_x_vals[i], heatmap_y_vals[j], height])
        end
    end

    # Normalize heatmap magnitude
    heatmap_mag .= heatmap_mag ./ maximum(heatmap_mag)

    # Extract balloon path from history
    balloon_x = history[1, :]  # Easting
    balloon_y = history[2, :]  # Northing
    balloon_time = LinRange(0, 1, size(history, 2))  # Time normalized to [0, 1]

    # Set up plot
    gr(legend=false, dpi=600)

    # Heatmap for wind magnitude
    heatmap(
        heatmap_x_vals, heatmap_y_vals, heatmap_mag,
        color=:viridis,
        title="Wind Field with Balloon Path at Height = $height m",
        xlabel="East [m]",
        ylabel="North [m]",
        colorbar=true
    )

    # Overlay quiver plot with scaled vectors
    quiver!(x_coords, y_coords, quiver=(u_scaled, v_scaled), color=:black)

    # Overlay balloon path color-coded by time
    scatter!(balloon_x, balloon_y, color=balloon_time, marker_z=balloon_time, ms=4, label="", c=:plasma)

    # Save and display plot
    Plots.savefig("./figures/WindTopDownWithPath.png")
    show()
    return nothing
end



# Runge-Kutta 4th-order integrator
function RK4_int(dt, state, fdot::Function)
    k1 = dt * fdot(state)
    k2 = dt * fdot(state + 0.5 * k1)
    k3 = dt * fdot(state + 0.5 * k2)
    k4 = dt * fdot(state + k3)
    return state + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
end