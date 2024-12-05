using PyPlot
# Example 3D trajectory
t = LinRange(0, 10, 100)
x = sin.(t)
y = cos.(t)
z = t

# 2D slices


# 3D trajectory
fig = figure(figsize=(10, 7))
ax = fig.add_subplot(111,projection="3d")
ax.plot(x, y, z, label="3D Trajectory", color="green")
xz_plane = ax.plot(x, z, label="xz-plane", color="blue")
xy_plane = ax.plot(x, y, label="xy-plane", color="red")
# Adding the 2D slices
ax.plot(x, zeros(length(x)), z, color="blue", linestyle="--", label="xz-slice")
ax.plot(x, y, zeros(length(x)), color="red", linestyle="--", label="xy-slice")

# Labels and Title
ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")
ax.set_title("3D Trajectory with xz and xy slices")

# Show the plot
legend()
show()