using LinearAlgebra, BenchmarkTools, PyPlot, Optimization, Plots
using Base.MathConstants: π

# Constants and parameters for the shaping functions
const D = 1000.0  # Characteristic downburst diameter
const zm = 0.025 * D  # Height of maximum horizontal velocity
const rm = 1.125 * D  # Radial distance of maximum horizontal velocity
const l = 10.5  # Scaling factor (to be calibrated)

# Shaping function parameters
const c1 = -0.16
const c2 = 1 / (1 + abs(c1))
const d = 2.1
const g = 0.5
const e = 0.3
const k = 0.1
const w = 1.0

# Shaping functions
function radial_shaping_function(r)
    if rm == 0
        error("Characteristic radius (rm) cannot be zero.")
    end
    r_normalized = r / rm
    return l * (r_normalized^d) * exp(-g * r_normalized^2) +
           e * exp(-k * r_normalized^2)
end
function vertical_shaping_function(z)
    if zm == 0
        error("Characteristic height (zm) cannot be zero.")
    end
    z_normalized = max(z / zm, 1e-6)  # Avoid zero or negative values
    return (z_normalized^(c2 - 1)) * exp(c1 * z_normalized^c2)
end

function vertical_velocity_function(r, z)
    if rm == 0 || zm == 0
        error("Characteristic dimensions (rm, zm) cannot be zero.")
    end
    r_normalized = r / rm
    z_normalized = max(z / zm, 1e-6)  # Avoid zero or negative values
    term1 = -l * g * exp(-g * r_normalized^2)
    term2 = c1 * c2 * exp(c1 * z_normalized^c2)
    return term1 * term2
end

function DownburstWindField(position::Vector{Float64})
    if length(position) != 3
        error("Position vector must have exactly 3 elements: [x, y, z].")
    end

    r = sqrt(position[1]^2 + position[2]^2)  # Radial distance
    z = position[3]  # Height

    u = radial_shaping_function(r) * vertical_shaping_function(z)
    v = 0.0  # Assume no tangential component for symmetry
    w = vertical_velocity_function(r, z)

    if isnan(u) || isnan(w)
        println("Warning: NaN detected in wind field components.")
        return [0.0, 0.0, 0.0]
    end

    return [u; v; w]
end


function WindField(pos)
    x,y,z = pos[1:3];
    W_N = 500*(z*z)  # North wind component
    W_E = 5    # East wind component
    W_U = 0   # Updraft component (f(x, y, z) if position-dependent)
    return [W_N; W_E; W_U]
end
using PyPlot

function plot_wind_direction_heatmap(domain_x, domain_z, resolution=50)
    # Define the grid over the domain
    x_vals = LinRange(domain_x[1], domain_x[2], resolution)
    z_vals = LinRange(domain_z[1], domain_z[2], resolution)

    # Create mesh grid
    x_coords = repeat(x_vals, 1, resolution)
    z_coords = repeat(z_vals', resolution, 1)

    # Compute wind field
    u = zeros(resolution, resolution)
    w = zeros(resolution, resolution)
    for i in 1:resolution
        for j in 1:resolution
            u[i, j], _, w[i, j] = DownburstWindField([x_coords[i, j], 0.0, z_coords[i, j]])
        end
    end

    # Plot using PyPlot
    figure(figsize=(12, 4))
    suptitle("Wind Field Heatmaps", fontsize=14)

    # First subplot for `u`
    subplot(1, 2, 1)
    imshow(u, extent=(domain_x[1], domain_x[2], domain_z[1], domain_z[2]), origin="lower", aspect="auto", cmap="viridis")
    colorbar()
    title("U Component")
    xlabel("x")
    ylabel("z")

    # Second subplot for `w`
    subplot(1, 2, 2)
    imshow(w, extent=(domain_x[1], domain_x[2], domain_z[1], domain_z[2]), origin="lower", aspect="auto", cmap="viridis")
    colorbar()
    title("W Component")
    xlabel("x")
    ylabel("z")

    # Save and show the plot
    show()
    PyPlot.savefig("./figures/WindHeatmap.png", dpi=300)
end


function plot_wind_field_3d(domain, resolution=10)
    x_vals = range(domain[1][1], domain[1][2], length=resolution)
    y_vals = range(domain[2][1], domain[2][2], length=resolution)
    z_vals = range(domain[3][1], domain[3][2], length=resolution)

    X, Y, Z = [], [], []
    U, V, W = [], [], []

    for x in x_vals
        for y in y_vals
            for z in z_vals
                u, v, w = DownburstWindField([x, y, z])
                push!(X, x)
                push!(Y, y)
                push!(Z, z)
                push!(U, u)
                push!(V, v)
                push!(W, w)
            end
        end
    end

    # Plot the 3D vector field
    Plots.quiver(
        X/1000, Y/1000, Z,
        quiver=(U, V, W),
        xlabel="x (km)", ylabel="y (km)", zlabel="z (m)", title="3D Wind Field",
        marker=:circle, linealpha=0.7
    )
    Plots.savefig("./figures/WindField2.png")
    Plots.quiver(
        X/1000,Z,
        quiver=(U,W);
        xlabel="x (km)", ylabel= "z (m)", title="2D Wind Field",
        marker=:circle, linealpha=0.7
    )
    show()
    Plots.savefig("./figures/WindField3.png")
    return nothing;
end

# Runge-Kutta 4th-order integrator
function RK4_int(dt, state, fdot::Function)
    k1 = dt * fdot(state)
    k2 = dt * fdot(state + 0.5 * k1)
    k3 = dt * fdot(state + 0.5 * k2)
    k4 = dt * fdot(state + k3)
    return state + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
end