using PyPlot
using PyCall

# Convert Julia types to numpy-compatible types
np = pyimport("numpy")  # Import numpy
X = np.linspace(1, 1, 100)
Y = np.linspace(1, 1, 100)
Z = np.ones((100, 100))
println(size(Z))

# Create 3D plot
fig = plt.figure()
ax = fig.add_subplot(projection="3d")

# Plot the 3D surface
# Plot projections of the contours for each dimension
ax.contourf(X, Y, Z, zdir="z", offset=-100, cmap="coolwarm")
ax.contourf(X, Y, Z, zdir="x", offset=-40, cmap="coolwarm")
ax.contourf(X, Y, Z, zdir="y", offset=40, cmap="coolwarm")

# Set axis limits and labels
ax.set(xlim=(-40, 40), ylim=(-40, 40), zlim=(-100, 100),
       xlabel="X", ylabel="Y", zlabel="Z")

plt.show()
