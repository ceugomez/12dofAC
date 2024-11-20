using LinearAlgebra, BenchmarkTools, PyPlot, Optimization
# Returns inertial wind vector as a function of inertial position (E, N, U)
function windField(pos)
    W_N = 10  # North wind component
    W_E = 0   # East wind component
    W_U = 1   # Updraft component (f(x, y, z) if position-dependent)
    return [W_N, W_E, W_U]
end

# Runge-Kutta 4th-order integrator
function RK4_int(dt, state, fdot)
    k1 = dt * fdot(state)
    k2 = dt * fdot(state + 0.5 * k1)
    k3 = dt * fdot(state + 0.5 * k2)
    k4 = dt * fdot(state + k3)
    return state + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
end