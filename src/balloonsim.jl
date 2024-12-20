include("generic.jl");
using Parameters, Atmosphere
# code to simulate pseudo-lagrangian drifter trajectories through a wind field
@with_kw mutable struct balloonParams
    # known from swenson et al; 
    # payload mass: 91g 
    # 125 L volume, ~ 0.125 m^3
    gravity::Float64 = 9.81; 
    ρa::Float64 = 1.225;            # kg⋅m³         (dynamically modified)
    ρh::Float64 = 0.166;            # kg⋅m³         (helium)
    Cd::Float64 = 0.5;              # drag coeff, dynamically modified
    r::Float64 = 0.310175245;       # m            
    vol::Float64 = 0.125;           # m³,          (swenson et.al.)
    m::Float64 = 0.092;             # kg           (swenson et.al)
    AGL::Float64 = 1611.7824;       # m, AGL->MSL conversion  (KBDU)
end
function balloonstatederiv(x::Vector{Float64})
    # Buoyancy force
        rho_z, mu, _ = atmospherefit(x[3]+prm.AGL)
        prm.ρa = rho_z
        fb = (prm.ρa - prm.ρh) * prm.gravity * prm.vol                  
        Fb = [0.0, 0.0, fb]  # Only in z direction

    # Drag force
        W = WindField(x[1:3])
        if (isnan(W[1]) || isnan(W[2]) || isnan(W[3]))
            println("error: NAN in wind field");
        end
        Vrel = x[4:6] - W  # Relative velocity
        # speed 
        rel_speed = norm(Vrel)
        prm.Cd = getDragCoeff(rel_speed,prm.ρa,mu);
        area = π * prm.r^2  
        Fd_mag = -prm.Cd * 0.5 * prm.ρa * rel_speed^2 * area
        Fd = Vrel./rel_speed * Fd_mag   # normlized direction*magnitude
    # Gravity force
        Fg = [0.0, 0.0, -prm.m * prm.gravity]
    # Acceleration: a = ∑F / m
        a = (Fb + Fd + Fg)./prm.m
    # State derivative: [velocity; acceleration]
    xdot = vcat(x[4:6], a)
    return xdot
end
function getDragCoeff(vrel::Float64, rho::Float64, mu::Float64)
    # calculate reynolds number
    RE = rho*abs(vrel)*2*prm.r/mu;
    Cd = 0.0;
    # calculate drag coefficient based on spherical drag model
    if (RE≤1)
        Cd = 24/RE
    end
    if (RE > 1 && RE ≤ 400)
        Cd = 24/(RE^0.646)
    end
    if (RE>400 && RE < 3*10^5)
        Cd = 0.5 
    end
    if (RE>3*10^5)
        Cd = 3.66*10^(-4)*(RE^0.4275)
    end
    return Cd
end

# EU quiver wind plotfunction plot_wind_direction_VF(domain_x, domain_z, resolution=50, arrow_scale=0.2)
# function plot_wind_direction_VF(domain_x, domain_z, history, resolution=50, arrow_scale=0.2)

#     # Generate x and z vectors
#     x_vals = LinRange(domain_x[1], domain_x[2], resolution)
#     z_vals = LinRange(domain_z[1], domain_z[2], resolution)

#     # Create meshgrid for x and z
#     x_coords = repeat(x_vals, 1, resolution)  # x repeated along rows
#     z_coords = repeat(z_vals', resolution, 1)  # z repeated along columns

#     # Initialize wind field components
#     u = zeros(resolution, resolution)
#     v = zeros(resolution, resolution)
#     w = zeros(resolution, resolution)

#     # Compute wind field values at each grid point
#     for i in 1:resolution
#         for j in 1:resolution
#             u[i, j], v[i,j], w[i, j] = WindField([x_coords[i, j], 0.0, z_coords[i, j]])
#         end
#     end

#     # Calculate magnitude of the wind field
#     mag = sqrt.(u .^ 2 + v .^ 2 + w .^ 2)

#     # Normalize vectors for quiver plot and scale them
#     u_scaled = arrow_scale .* (u ./ (mag .+ eps()))
#     w_scaled = arrow_scale .* (w ./ (mag .+ eps()))

#     # Create the plot
#     fig, ax = subplots()
#     im = ax.imshow(
#         mag',
#         extent=(domain_x[1], domain_x[2], domain_z[1], domain_z[2]),
#         origin="lower",
#         cmap="viridis",
#         aspect="auto"
#     )
#     colorbar(im, ax=ax, label="Wind Magnitude [m/s]")
#     downscale = 10
#     ax.quiver(
#         x_coords[1:downscale:end,1:downscale:end], z_coords[1:downscale:end,1:downscale:end], u_scaled[1:downscale:end,1:downscale:end], w_scaled[1:downscale:end,1:downscale:end],
#         color="black", scale=1/arrow_scale
#     )
#     ax.set_title("Wind Field (Vertical Slice at N=0)")
#     ax.set_xlabel("East [m]")
#     ax.set_ylabel("Up [m]")
#     balloon_x = history[1, :]  # Easting
#     balloon_y = history[2, :]  # Northing
#     balloon_z = history[3, :]  # Altitude
#     ax.plot(
#             balloon_x, balloon_z, 
#             label="Trajectory", 
#             color="Black", 
#             linewidth=1, 
#             zorder=8  
#         )
#         ax.scatter(
#             balloon_x[1], balloon_z[1],
#             s=[75],
#             alpha=1.0,
#             label="Release",
#             color="green",
#             zorder=10
#         )
#         ax.scatter(
#             balloon_x[end],balloon_z[end],
#             s=[75],
#             alpha=1.0,
#             label="Final",
#             color="Red",
#             zorder=10
#         )
#         legend()

#     # Save the plot
#     fig.savefig("./figures/WindVF_EU.png")
#     show()
#     return nothing
# end

# function plot_wind_top_down(domain_x, domain_y, height=50.0, resolution=50, arrow_scale=0.2)
#     # Generate x and y vectors
#     x_vals = LinRange(domain_x[1], domain_x[2], resolution)
#     y_vals = LinRange(domain_y[1], domain_y[2], resolution)
#     # Meshgrid
#     x_coords = repeat(x_vals, 1, resolution)
#     y_coords = repeat(y_vals', resolution, 1)

#     # Compute wind field at specified height
#     u = zeros(resolution, resolution)
#     v = zeros(resolution, resolution)
#     for i in 1:resolution
#         for j in 1:resolution
#             u[i, j], v[i, j], _ = WindField([x_coords[i, j], y_coords[i, j], height])
#         end
#     end

#     # Calculate wind magnitude
#     mag = sqrt.(u .^ 2 + v .^ 2)

#     # Normalize the vectors and scale them for the quiver plot
#     u_scaled = arrow_scale .* (u ./ (mag .+ eps()))  
#     v_scaled = arrow_scale .* (v ./ (mag .+ eps()))

#     # Set up plot
#     gr(legend=false, dpi=600)

#     # Heatmap for wind magnitude
#     heatmap(
#         x_vals, y_vals, mag,
#         color=:viridis,
#         title="Wind Field at Height = $height m",
#         xlabel="East [m]",
#         ylabel="North [m]"
#     )

#     # Overlay quiver plot with scaled vectors
#     quiver!(x_coords[1:5:end], y_coords[1:5:end], quiver=(u_scaled[1:5:end,1:5:end], v_scaled[1:5:end,1:5:end]), color=:black)

#     # Save and display plot
#     Plots.savefig("./figures/WindTopDown.png")
#     show()
#     return nothing
# end

# function plot_wind_top_down_with_path(domain_x, domain_y, history::Matrix{Float64}, height=50.0, resolution=50, arrow_scale=0.2)
#     # Generate x and y vectors
#     downscale = 20
#     x_vals = LinRange(domain_x[1], domain_x[2], resolution)
#     y_vals = LinRange(domain_y[1], domain_y[2], resolution)
#     # Meshgrid
#     x_coords = repeat(x_vals, 1, resolution)
#     y_coords = repeat(y_vals', resolution, 1)

#     # Compute wind field at specified height
#     u = zeros(resolution, resolution)
#     v = zeros(resolution, resolution)
#     w = zeros(resolution, resolution)
#     for i in 1:resolution
#         for j in 1:resolution
#             u[i, j], v[i, j], w[i,j] = WindField([x_coords[i, j], y_coords[i, j], height])
#         end
#     end

#     # Calculate wind magnitude
#     mag = sqrt.(u .^ 2 + v .^ 2 + w .^ 2)

#     # Normalize the vectors and scale them for the quiver plot
#     u_scaled = arrow_scale .*u 
#     v_scaled = arrow_scale .*v
#     # Extract balloon path from history
#     balloon_x = history[1, :]  # Easting
#     balloon_y = history[2, :]  # Northing

#     # Create the plot
#     fig = figure(figsize=(12, 9))
#     ax = fig.add_subplot(111)
#     im = ax.imshow(
#         mag',
#         extent=(domain_x[1], domain_x[2], domain_y[1], domain_y[2]),
#         origin="lower",
#         cmap="viridis",
#         aspect="auto"
#     )
#     colorbar(im, ax=ax, label="Wind Magnitude [m/s]")
#     downscale = 10
#     ax.quiver(
#         x_coords[1:downscale:end,1:downscale:end], y_coords[1:downscale:end,1:downscale:end], u_scaled[1:downscale:end,1:downscale:end], v_scaled[1:downscale:end,1:downscale:end],
#         color="black", scale=1/arrow_scale
#     )
#     ax.set_title("Wind Field (Vertical Slice at $height m)")
#     ax.set_xlabel("East [m]")
#     ax.set_ylabel("North [m]")
#     ax.plot(balloon_x,balloon_y, color="black", label="Trajectory")
#     ax.set_xlim(domain_x[1:2])
#     ax.set_ylim(domain_y[1:2])
#     ax.scatter(
#             balloon_x[1], balloon_y[1],
#             s=[75],
#             alpha=1.0,
#             label="Release",
#             color="green",
#             zorder=10
#         )
#     legend()

#     # Save the plot
#     fig.savefig("./figures/WindVF_EN.png")
#     show()
#     return nothing
# end

# function PlotWindField3d(domain_x, domain_y, domain_z, history, height=50, resolution=50, arrow_scale=0.2)
#     np = pyimport("numpy")
#     downscale = 10;
#     arrow_scale = 100
#     # Generate vectors for x, y, and z
#         x_vals = LinRange(domain_x[1], domain_x[2], resolution)
#         y_vals = LinRange(domain_y[1], domain_y[2], resolution)
#         z_vals = LinRange(domain_z[1], domain_z[2], resolution)

#     # xy plane @ z=height
#         uxy = zeros(resolution, resolution)
#         vxy = zeros(resolution, resolution)
#         for i in 1:resolution
#             for j in 1:resolution
#                 uxy[i, j], vxy[i, j], _ = WindField([x_vals[i], y_vals[j], height])
#             end
#         end
#         xymag = sqrt.(uxy.^2 .+ vxy.^2)
#     # xz plane @ y=0
#         uxz = zeros(resolution, resolution)
#         wxz = zeros(resolution, resolution)
#         for i in 1:resolution
#             for j in 1:resolution
#                 uxz[i, j], _, wxz[i, j] = WindField([x_vals[i], 0.0, z_vals[j]])
#             end
#         end
#         xzmag = sqrt.(uxz.^2 .+ wxz.^2)
#     # yz plane @ x=0
#         vyz = zeros(resolution, resolution)
#         wyz = zeros(resolution, resolution)
#         for i in 1:resolution
#             for j in 1:resolution
#                 _, vyz[i, j], wyz[i, j] = WindField([0.0, y_vals[i], z_vals[j]])
#             end
#         end
#         yzmag = sqrt.(vyz.^2 .+ wyz.^2)    
#     # Create the figure
#         fig = figure(figsize=(12, 9))
#         ax = fig.add_subplot(111, projection="3d")
#         ax.grid("False")
#     # Labels and title
#         plt.grid("None")
#         ax.set_xlabel("Easting")
#         ax.set_ylabel("Northing")
#         ax.set_zlabel("Altitude")
#         ax.set_xlim(domain_x[1:2])
#         ax.set_ylim(domain_y[1:2])
#         ax.set_zlim(domain_z[1:2])
#         ax.set_title("Downburst Wind Field")
#     # Contours   
#         # Plot xy plane contour
#                 X, Y = np.meshgrid(x_vals, y_vals)  # Grid for xy plane
#                 ax.contourf(
#                     X, Y, xymag,
#                     zdir="z",    # Projection direction
#                     offset=domain_z[1],  # Bottom of z-axis
#                     cmap="viridis",
#                     alpha=0.6,
#                     zorder=1
#                 )
#         # Plot xz plane contour plot
#                 X, Z = np.meshgrid(x_vals, z_vals)  # Grid for xz plane
#                 ax.contourf(
#                     X, xzmag', Z,
#                     zdir="y",  
#                     cmap="viridis",
#                     alpha=0.6,
#                     zorder=1,
#                     offset=domain_y[2]  # Place on y=center_y
#                 )
#         # Plot yz plane contour
#                 Y, Z = np.meshgrid(y_vals, z_vals)  # Grid for xz plane
#                 ax.contourf(
#                     yzmag', Y, Z,
#                     zdir="x",  
#                     cmap="viridis",
#                     alpha=0.6,
#                     zorder=1,
#                     offset=domain_x[1]  # Place on x=center_x
#                 ) 
#     # Quivers
#         # --- Quiver for xy plane ---
#         # Downscale quiver
#         x_quiver = x_vals[1:downscale:end]
#         y_quiver = y_vals[1:downscale:end]
#         uxy_downscaled = uxy[1:downscale:end, 1:downscale:end]
#         vxy_downscaled = vxy[1:downscale:end, 1:downscale:end]

#         # Crop the edges by excluding the outermost vectors
#         x_quiver = x_quiver[2:end-1]  # Remove outermost x values
#         y_quiver = y_quiver[2:end-1]  # Remove outermost y values
#         uxy_downscaled = uxy_downscaled[2:end-1, 2:end-1]  # Remove outermost vectors
#         vxy_downscaled = vxy_downscaled[2:end-1, 2:end-1]  # Remove outermost vectors

#         # Create grid for xy plane quiver
#         x_quiver_plane, y_quiver_plane = repeat(x_quiver, 1, length(y_quiver)), repeat(y_quiver', length(x_quiver), 1)
#         z_quiver_plane = 2 * ones(size(x_quiver_plane))

#         ax.quiver(
#             x_quiver_plane, y_quiver_plane, z_quiver_plane,
#             uxy_downscaled, vxy_downscaled, zeros(size(uxy_downscaled)),
#             color="black", length=arrow_scale, normalize="False", zorder=10
#         )

#         # --- Quiver for xz plane ---
#         # Downscale quiver for xz plane
#         x_quiver = x_vals[1:downscale:end]
#         z_quiver = z_vals[1:downscale:end]
#         uxz_downscaled = uxz[1:downscale:end, 1:downscale:end]
#         wxz_downscaled = wxz[1:downscale:end, 1:downscale:end]

#         # Crop the edges
#         x_quiver = x_quiver[2:end-1]
#         z_quiver = z_quiver[2:end-1]
#         uxz_downscaled = uxz_downscaled[2:end-1, 2:end-1]
#         wxz_downscaled = wxz_downscaled[2:end-1, 2:end-1]

#         # Create grid for xz plane quiver
#         x_quiver_plane, z_quiver_plane = repeat(x_quiver, 1, length(z_quiver)), repeat(z_quiver', length(x_quiver), 1)
#         y_quiver_plane = 997.0 * ones(size(x_quiver_plane))

#         ax.quiver(
#             x_quiver_plane, y_quiver_plane, z_quiver_plane,
#             uxz_downscaled, zeros(size(uxz_downscaled)), wxz_downscaled,
#             color="black", length=arrow_scale, normalize="False", zorder=10
#         )

#         # --- Quiver for yz plane ---
#         # Downscale quiver for yz plane
#         y_quiver = y_vals[1:downscale:end]
#         z_quiver = z_vals[1:downscale:end]
#         vyz_downscaled = vyz[1:downscale:end, 1:downscale:end]
#         wyz_downscaled = wyz[1:downscale:end, 1:downscale:end]

#         # Crop the edges
#         y_quiver = y_quiver[2:end-1]
#         z_quiver = z_quiver[2:end-1]
#         vyz_downscaled = vyz_downscaled[2:end-1, 2:end-1]
#         wyz_downscaled = wyz_downscaled[2:end-1, 2:end-1]

#         # Create grid for yz plane quiver
#         y_quiver_plane, z_quiver_plane = repeat(y_quiver, 1, length(z_quiver)), repeat(z_quiver', length(y_quiver), 1)
#         x_quiver_plane = -997 * ones(size(y_quiver_plane))

#         ax.quiver(
#             x_quiver_plane, y_quiver_plane, z_quiver_plane,
#             zeros(size(vyz_downscaled)), vyz_downscaled, wyz_downscaled,
#             color="black", length=arrow_scale, normalize="False", zorder=10
#         )
#     # Plot balloon path from history
#         balloon_x = history[1, :]  # Easting
#         balloon_y = history[2, :]  # Northing
#         balloon_z = history[3, :]  # Altitude
#         ax.plot(
#             balloon_x, balloon_y, balloon_z, 
#             label="3D Trajectory", 
#             color="Black", 
#             linewidth=5, 
#             zorder=8  
#         )
#         ax.scatter(
#             balloon_x[1], balloon_y[1], balloon_z[1],
#             s=[75],
#             alpha=1.0,
#             label="Release",
#             color="green",
#             zorder=10
#         )
#         ax.scatter(
#             balloon_x[end], balloon_y[end], balloon_z[end],
#             s=[75],
#             alpha=1.0,
#             label="Final",
#             color="Red",
#             zorder=10
#         )
#         legend()
#         plt.savefig("./figures/3dplot.png", dpi=300)
#         plt.show()

#     return nothing
# end

# function visualize_field_xyz(history::Matrix{Float64})
#     # Validate input dimensions
#     if size(history, 1) < 3
#         error("Input matrix must have at least 3 rows representing E, N, and Up.")
#     end

#     # Create subplots for the tracks
#     plot1 = Plots.plot(history[1, :]/1000, history[2, :]/1000, xlabel="Easting [m]", ylabel="Northing [m]", title="NE Track")
#     plot2 = Plots.plot(history[1, :]/1000, history[3, :], xlabel="Easting [m]", ylabel="Up [m]", title="EU Track")
#     plot3 = Plots.plot(history[2, :]/1000, history[3, :], xlabel="Northing [m]", ylabel="Up [m]", title="3D Track")

#     # Combine the plots into a grid layout
#     layout = @layout [a b; c]
#     Plots.plot(plot1, plot2, plot3, layout=layout, size=(900, 600))

#     # Save the figure
#     Plots.savefig("./figures/Track.png")
#     show()
#     return nothing
# end

# function visualize_field_t(history::Matrix{Float64}, dt::Float64, tf::Float64)
#     # Validate input dimensions
#     @assert size(history, 1) == 3 "Input 'history' must have 3 rows (Easting, Northing, Up)."

#     # Create time vector in minutes
#     t = LinRange(0, size(history, 2) / (60/dt), size(history, 2))

#     # Set up the figure
#     plt = figure(figsize=(12, 4))
#     plt.suptitle("Balloon Time History", fontsize=14)

#     # Plot Easting vs. Time
#     ax1 = plt.add_subplot(131)
#     ax1.plot(t, history[1, :]/1000, label="Easting", color="blue")
#     ax1.set_xlabel("Time [min]")
#     ax1.set_ylabel("Easting [km]")
#     ax1.grid(true)
#     ax1.legend()

#     # Plot Northing vs. Time
#     ax2 = plt.add_subplot(132)
#     ax2.plot(t, history[2, :]/1000, label="Northing", color="green")
#     ax2.set_xlabel("Time [min]")
#     ax2.set_ylabel("Northing [km]")
#     ax2.grid(true)
#     ax2.legend()

#     # Plot Up vs. Time
#     ax3 = plt.add_subplot(133)
#     ax3.plot(t, history[3, :], label="Up", color="red")
#     ax3.set_xlabel("Time [min]")
#     ax3.set_ylabel("Up [m]")
#     ax3.grid(true)
#     ax3.legend()

#     # Adjust layout and save the figure
#     plt.tight_layout(rect=[0, 0, 1, 0.95])  # Leave space for the title
#     PyPlot.savefig("./figures/TTrack.png", dpi=300)
#     return nothing
# end
# # runtime loop
# begin
#     # Initialize parameter structure
#     global prm = balloonParams()

#     # Simulation parameters
#     tstep = 0.01  # Time step in seconds
#     tf = 1000.0 #5500.0  # Total simulation time in seconds
#     noIdx = Int64(tf / tstep) + 1  # Number of time steps including initial state

#     # Initialize history array (6 state variables: x, y, z, vx, vy, vz)
#     history = zeros(Float64, 6, noIdx)

#     # Set initial state 
#     history[:, 1] = [10.0, 25.0, 70.0, 0.0, 0.0, 0.0]  # [x, y, z, vx, vy, vz]

#     # Simulation loop
#     for i in 1:(noIdx - 1)
#         # Get current state
#         stateK = history[:, i]

#         # Integrate to find next state
#         stateKp1 = RK4_int(tstep, stateK, balloonstatederiv)

#         # Store next state in history
#         history[:, i + 1] = stateKp1
#     end

#     println("Simulation complete.")

#     # Define the domain and resolution
#     domain_x = (-1000, 1000)  # x domain in meters
#     domain_y = (-1000, 1000)  # y domain in meters
#     domain_z = (0,1500)     # z domain in meters
#     resolution = 100         # Number of grid points
#     height_z = 40.0;

#     # Plot the wind direction
#     #plotFancy(domain_x,domain_y,domain_z,history,50.0,resolution,50.0)
#     #PleasePleasePlease(domain_x,domain_y,domain_z,history,50.0,resolution,50.0)
#     plot_wind_direction_VF(domain_x, domain_z, history, resolution, 0.2)
#     #plot_wind_top_down(domain_x, domain_y, height_z, resolution, 0.2)
#     plot_wind_top_down_with_path(domain_x, domain_y, history, height_z, resolution, 0.075)
#     visualize_field_t(history[1:3, :],tstep,tf)  
#     #visualize_field_xyz(history[1:3, :])
#     println("Plotted")
# end
