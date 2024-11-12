# cgf cego6160@colorado.edu
# 12DOF aircraft simulation in Julia
struct parameters
    

end



#state derivatives as a function of state, wind, control, and aircraft parameters
function ACstatederiv(state,w,u, param)
    # state: [N E D u v w phi theta psi p q r] 
    inertial_position_earth = state[1:3];
    inertial_velocity_body = state[4:6];
    body_angles = state[7:9];
    body_rates = state[10:12];






    return statedot;
end
# returns the state derivative of a superpressure balloon in a wind field
function balloonstatederiv(state)
    return state;
end
# returns inertial wind vector as a function of inertial position (E,N,U)
function windField(pos)
    W_N = 10; # north wind
    W_E = 0;  # East wind
    W_U = 1;  # updraft component f(xyz)
end
# Runge-Kutta 4th-order integrator
function RK4_int(dt, state, fdot)

end
# Construct rotation matrix 
function rotate(c1,c2,angles)
begin
    println("Starting Sim");
    # state: [N E D u v w phi theta psi p q r] 


end