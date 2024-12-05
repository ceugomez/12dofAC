using LinearAlgebra, FFMPEG, BenchmarkTools, PyPlot, Optimization, Plots
using Base.MathConstants: π
# Constants and parameters for the shaping functions
const D = 2500.0       # Characteristic downburst diameter [m]

# params from Vicroy NASA TM 104053
const rp = 650.0; # Radius of peak horizontal wind [m]
const zm = 60;  # altitude of peak vertical wind [m]
const um = 10;  # magnitude of maximum outflow [ms⁻¹]
const α = 2;    # const
const c1 = -0.22
const c2 = -2.75
const λ = (2*um)/(rp*(exp(c1)-exp(c2))*exp(1/(2*α)));   # Vicroy, eq. A1
const β = (2*rp^(2*α))^(1/(2*α))

function WindField(pos::Vector{Float64})
    x, y, z = pos[1:3]
    

    u = (λ*x/2)*
        (exp(c1*(z/zm))-exp(c2*(z/zm)))*
        exp(
            (2-((x^2+y^2)^α)/(rp^(2α)))
            /(2*α)
            )
    v = (λ*y/2)*
        (exp(c1*(z/zm))-exp(c2*(z/zm)))*
        exp(
            (2-((x^2+y^2)^α)/(rp^(2α)))
            /(2*α)
            )
    w = -λ*(
            (zm/c1)*(exp(c1*(z/zm))-1)-(zm/c2)*exp(c2*(z/zm)-1)
        )*(
            1-((x^2+y^2)^α)/((2*rp)^(2*α))
        )*exp(
            (2-((x^2+y^2)^α)/(rp^(2α)))
            /(2*α)
        )
    return [u,v,w]
    # # Horizontal wind component (u)
    # u = (λ * r / 2) *
    #     (
    #         exp(ζ1 * (z / zm)) - exp(ζ2 * (z / zm))
    #     ) *
    #     exp(
    #         (1 - (r^2 / β^2)^α) / α
    #     )

    # # Vertical wind component (w)
    # w = (
    #         (zm / ζ1) * (exp(ζ1 * (z / zm)) - 1)
    #         -
    #         (zm / ζ2) * (exp(ζ2 * (z / zm)) - 1)
    #     ) *
    #     (
    #         1 - (r^2 / β^2)^α
    #     ) *
    #     exp(
    #         (1 - (r^2 / β^2)^α) / α
    #     )

    # # Cylindrical to Cartesian conversion
    # if E != 0 || N != 0
    #     θ = atan(N, E)
    # else
    #     θ = π / 4
    # end

    # return [u * cos(θ), u * sin(θ), w]
end

# generic wind field 
function GenericWindField(pos)
    E, N, U = pos[1:3]
    W_N = U / 10
    W_E = U / 10
    W_U = 0
    return [W_E; W_N; W_U]
end
# EU quiver wind plotfunction plot_wind_direction_VF(domain_x, domain_z, resolution=50, arrow_scale=0.2)
function plot_wind_direction_VF(domain_x, domain_z, resolution=50, arrow_scale=0.2)

    # Generate x and z vectors
    x_vals = LinRange(domain_x[1], domain_x[2], resolution)
    z_vals = LinRange(domain_z[1], domain_z[2], resolution)

    # Create meshgrid for x and z
    x_coords = repeat(x_vals, 1, resolution)  # x repeated along rows
    z_coords = repeat(z_vals', resolution, 1)  # z repeated along columns

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
    u_scaled = arrow_scale .* (u ./ (mag .+ eps()))
    w_scaled = arrow_scale .* (w ./ (mag .+ eps()))

    # Create the plot
    fig, ax = subplots()
    im = ax.imshow(
        mag',
        extent=(domain_x[1], domain_x[2], domain_z[1], domain_z[2]),
        origin="lower",
        cmap="viridis",
        aspect="auto"
    )
    colorbar(im, ax=ax, label="Wind Magnitude [m/s]")
    downscale = 10
    ax.quiver(
        x_coords[1:downscale:end,1:downscale:end], z_coords[1:downscale:end,1:downscale:end], u_scaled[1:downscale:end,1:downscale:end], w_scaled[1:downscale:end,1:downscale:end],
        color="black", scale=1/arrow_scale
    )
    ax.set_title("Wind Field (Vertical Slice at N=0)")
    ax.set_xlabel("East [m]")
    ax.set_ylabel("Up [m]")

    # Save the plot
    fig.savefig("./figures/WindVF_EU.png")
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
    quiver!(x_coords[1:5:end], y_coords[1:5:end], quiver=(u_scaled[1:5:end,1:5:end], v_scaled[1:5:end,1:5:end]), color=:black)

    # Save and display plot
    Plots.savefig("./figures/WindTopDown.png")
    show()
    return nothing
end

function plot_wind_top_down_with_path(domain_x, domain_y, history::Matrix{Float64}, height=50.0, resolution=50, arrow_scale=0.2)
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
    u_scaled = arrow_scale .*u 
    v_scaled = arrow_scale .*v
    # Extract balloon path from history
    balloon_x = history[1, :]  # Easting
    balloon_y = history[2, :]  # Northing

    # Create the plot
    fig, ax = subplots()
    im = ax.imshow(
        mag',
        extent=(domain_x[1], domain_x[2], domain_y[1], domain_y[2]),
        origin="lower",
        cmap="viridis",
        aspect="auto"
    )
    colorbar(im, ax=ax, label="Wind Magnitude [m/s]")
    downscale = 10
    ax.quiver(
        x_coords[1:downscale:end,1:downscale:end], y_coords[1:downscale:end,1:downscale:end], u_scaled[1:downscale:end,1:downscale:end], v_scaled[1:downscale:end,1:downscale:end],
        color="black", scale=1/arrow_scale
    )
    ax.set_title("Wind Field (Vertical Slice at $height m)")
    ax.set_xlabel("East [m]")
    ax.set_ylabel("Up [m]")
    ax.plot(balloon_x,balloon_y)
    ax.legend("Balloon Path")

    # Save the plot
    fig.savefig("./figures/WindVF_EN.png")
    show()
    return nothing
end
function plotFancy(domain_x, domain_y, domain_z, history, height=50, resolution=50, arrow_scale=0.2)
    downscale = 10;
    # Generate vectors for x, y, and z
    x_vals = LinRange(domain_x[1], domain_x[2], resolution)
    y_vals = LinRange(domain_y[1], domain_y[2], resolution)
    z_vals = LinRange(domain_z[1], domain_z[2], resolution)

    # Compute fields for yz, xz, and xy planes
    # yz plane @ x=-1000
    vyz = zeros(resolution, resolution)
    wyz = zeros(resolution, resolution)
    for i in 1:resolution
        for j in 1:resolution
            _, vyz[i, j], wyz[i, j] = WindField([0.0, y_vals[i], z_vals[j]])
        end
    end
    yzmag = sqrt.(vyz.^2 .+ wyz.^2)

    # xz plane @ y=-1000
    uxz = zeros(resolution, resolution)
    wxz = zeros(resolution, resolution)
    for i in 1:resolution
        for j in 1:resolution
            uxz[i, j], _, wxz[i, j] = WindField([x_vals[i], 0.0, z_vals[j]])
        end
    end
    xzmag = sqrt.(uxz.^2 .+ wxz.^2)

    # xy plane @ z=height
    uxy = zeros(resolution, resolution)
    vxy = zeros(resolution, resolution)
    for i in 1:resolution
        for j in 1:resolution
            uxy[i, j], vxy[i, j], _ = WindField([x_vals[i], y_vals[j], height])
        end
    end
    xymag = sqrt.(uxy.^2 .+ vxy.^2)
    # Extract balloon path from history
    balloon_x = history[1, :]  # Easting
    balloon_y = history[2, :]  # Northing
    balloon_z = history[3, :]  # Altitude
    # Create the figure
        fig = figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection="3d")
        ax[:set_axis_off]()
        # Labels and title
        ax.grid("False")
        plt.grid("None")
        ax.set_xlabel("Easting")
        ax.set_ylabel("Northing")
        ax.set_zlabel("Altitude")
        ax.set_xlim(domain_x[1:2])
        ax.set_ylim(domain_y[1:2])
        ax.set_zlim(domain_z[1:2])
        ax.set_title("Downburst Wind Field")
    # Plot balloon path
        ax.plot(balloon_x, balloon_y, balloon_z, label="3D Trajectory", color="green")

    # Plot yz plane @ x=-1000
        y_plane, z_plane = repeat(y_vals, 1, resolution), repeat(z_vals', resolution, 1)
        x_plane = -999 * ones(size(y_plane))
        ax.plot_surface(
            x_plane, y_plane, z_plane,
            rstride=1, cstride=1, facecolors=PyPlot.cm.viridis(yzmag./maximum(xzmag)),
            shade="True", alpha=0.4
        )
    # Plot xz plane @ y=-1000
        x_plane, z_plane = repeat(x_vals, 1, resolution), repeat(z_vals', resolution, 1)
        y_plane = 999.0 * ones(size(x_plane))
        ax.plot_surface(
            x_plane, y_plane, z_plane,
            rstride=1, cstride=1, facecolors=PyPlot.cm.viridis(xzmag./maximum(xzmag)),
            shade="True", alpha=0.4
        )
    # Plot xy plane @ z=1
        x_plane, y_plane = repeat(x_vals, 1, resolution), repeat(y_vals', resolution, 1)
        z_plane = zeros(size(x_plane))
        ax.plot_surface(
            x_plane, y_plane, z_plane,
            rstride=1, cstride=1, facecolors=PyPlot.cm.viridis(xymag./maximum(xymag)),
            shade="True", alpha=0.4
        )

    # Quiver for xy plane
        x_quiver = x_vals[1:downscale:end]
        y_quiver = y_vals[1:downscale:end]
        uxy_downscaled = uxy[1:downscale:end, 1:downscale:end]
        vxy_downscaled = vxy[1:downscale:end, 1:downscale:end]
        x_quiver_plane, y_quiver_plane = repeat(x_quiver, 1, length(y_quiver)), repeat(y_quiver', length(x_quiver), 1)
        z_quiver_plane = 2*ones(size(x_quiver_plane))
        ax.quiver(
            x_quiver_plane, y_quiver_plane, z_quiver_plane,
            uxy_downscaled, vxy_downscaled, zeros(size(uxy_downscaled)),
            color="black", length=arrow_scale, normalize="False"
        )
    # Quiver for xz plane
        x_quiver = x_vals[1:downscale:end]
        z_quiver = z_vals[1:downscale:end]
        uxz_downscaled = uxz[1:downscale:end, 1:downscale:end]
        wxz_downscaled = wxz[1:downscale:end, 1:downscale:end]
        x_quiver_plane, z_quiver_plane = repeat(x_quiver, 1, length(z_quiver)), repeat(z_quiver', length(x_quiver), 1)
        y_quiver_plane = 997.0 * ones(size(x_quiver_plane))
        ax.quiver(
            x_quiver_plane, y_quiver_plane, z_quiver_plane,
            uxz_downscaled, zeros(size(uxz_downscaled)), wxz_downscaled,
            color="black", length=arrow_scale, normalize="False"
        )
    # Quiver for yz plane

        y_quiver = y_vals[1:downscale:end]
        z_quiver = z_vals[1:downscale:end]
        vyz_downscaled = vyz[1:downscale:end, 1:downscale:end]
        wyz_downscaled = wyz[1:downscale:end, 1:downscale:end]
        y_quiver_plane, z_quiver_plane = repeat(y_quiver, 1, length(z_quiver)), repeat(z_quiver', length(y_quiver), 1)
        x_quiver_plane = -997 * ones(size(y_quiver_plane))
        ax.quiver(
            x_quiver_plane, y_quiver_plane, z_quiver_plane,
            zeros(size(vyz_downscaled)), vyz_downscaled, wyz_downscaled,
            color="black", length=arrow_scale, normalize="False"
        )
    PyPlot.grid("False")
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