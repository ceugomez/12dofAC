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
begin
    println("Starting Simulation")
end
