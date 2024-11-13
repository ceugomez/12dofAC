# cgf cego6160@colorado.edu
# 12DOF aircraft simulation in Julia

# Define parameters struct (currently empty; add fields as needed)
struct Parameters
    # Define any necessary fields for aircraft parameters here
end

# Define aircraft state structure
struct AcState
    pe::Vector{Vector{Float64}}    # Position in Earth coordinates
    eang::Vector{Vector{Float64}}  # Euler angles
    vbe::Vector{Vector{Float64}}   # Velocity in the body frame
    veang::Vector{Vector{Float64}} # Angular velocity in the body frame
end

# State derivatives as a function of state, wind, control, and aircraft parameters
function ACstatederiv(state::AcState, w, u, param::Parameters)
    # state: [N E D u v w phi theta psi p q r]
    inertial_position_earth = state.pe
    inertial_velocity_body = state.vbe   # Fixed typo: should match field `vbe`
    body_angles = state.eang
    body_rates = state.veang
    
    # Placeholder for state derivatives computation (to be implemented)
    statedot = nothing # Change `nothing` to the actual derivative computation
    
    return statedot
end

# Function to compute the rotation matrix from Euler angles
function Rbe(eang::Vector{Float64})
    Rbe = [cos(eang[2])*cos(eang[1]) -sin(eang[1]) 0;  #Not done
           sin(eang[2])*cos(eang[1]) cos(eang[1])  0;
           0                      0                1]
    return Rbe
end

# State derivative of a balloon in a wind field (to be implemented)
function balloonstatederiv(state)
    # Placeholder for balloon state derivative computation
    return state
end

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
begin
    println("Starting Simulation")
end
