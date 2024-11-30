using LinearAlgebra, BenchmarkTools, PyPlot, Optimization
using Base.MathConstants: π

# Constants and parameters for the shaping functions
const D = 1000.0  # Characteristic downburst diameter
const zm = 0.025 * D  # Height of maximum horizontal velocity
const rm = 1.125 * D  # Radial distance of maximum horizontal velocity
const l = 0.5  # Scaling factor (to be calibrated)

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
    r_normalized = r / rm
    return l * (r_normalized^d) * exp(-g * r_normalized^2) +
           e * exp(-k * r_normalized^2)
end

function vertical_shaping_function(z)
    z_normalized = z / zm
    return (z_normalized^(c2 - 1)) * exp(c1 * z_normalized^c2)
end

function vertical_velocity_function(r, z)
    r_normalized = r / rm
    z_normalized = z / zm
    term1 = -l * g * exp(-g * r_normalized^2)
    term2 = c1 * c2 * exp(c1 * z_normalized^c2)
    return term1 * term2
end

# Wind field calculation
function DownburstWindField(position::Vector{Float64})
    r = sqrt(position[1]^2 + position[2]^2)  # Radial distance
    println(r)
    z = position[3]  # Height
    println(z)
    u = radial_shaping_function(r) * vertical_shaping_function(z)
    println(u)
    sleep(5)
    v = 0.0  # Assume no tangential component for symmetry
    w = vertical_velocity_function(r, z)
    if (isnan(u))
        println("u isNAN!")
        return [0,0,0]
    end
    return [u, v, w]
end

function WindField(pos)
    x,y,z = pos[1:3];
    W_N = 5*(z*z)  # North wind component
    W_E = 0    # East wind component
    W_U = 0   # Updraft component (f(x, y, z) if position-dependent)
    return [W_N; W_E; W_U]
end

using Plots

# Visualization of the wind field
function plot_wind_field_2d(domain_x, domain_z, resolution=20)
    x_vals = range(domain_x[1], domain_x[2], length=resolution)
    z_vals = range(domain_z[1], domain_z[2], length=resolution)
    u_vals = zeros(resolution, resolution)
    w_vals = zeros(resolution, resolution)

    for (i, x) in enumerate(x_vals)
        for (j, z) in enumerate(z_vals)
            u, _, w = WindField([x, 0.0, z])  # 2D assumes y=0
            u_vals[i, j] = u
            w_vals[i, j] = w
        end
    end

    # Plot the vector field
    quiver(
        x_vals, z_vals,
        quiver=(u_vals', w_vals'),
        xlabel="x (m)", ylabel="z (m)", title="2D Wind Field",
        aspect_ratio=:equal, legend=false
    )
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
                u, v, w = WindField([x, y, z])
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
    quiver(
        X, Y, Z,
        quiver=(U, V, W),
        xlabel="x (m)", ylabel="y (m)", zlabel="z (m)", title="3D Wind Field",
        marker=:circle, linealpha=0.7
    )
end

# Runge-Kutta 4th-order integrator
function RK4_int(dt, state, fdot::Function)
    k1 = dt * fdot(state)
    k2 = dt * fdot(state + 0.5 * k1)
    k3 = dt * fdot(state + 0.5 * k2)
    k4 = dt * fdot(state + k3)
    return state + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
end