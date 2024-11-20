include("generic.jl");
# State derivative of a balloon in a wind field (to be implemented)
mutable struct balloonParams
    g = 9.81
    ρa::Float64 = 1.225 # kg/m^3        !! Modifiable based on state
    ρh::Float64 = 0.166 # kg/m^3, helium
    Cd::Float64 = 0.01  # drag coeff    !!Placeholder
    r::Float64 = 1;     # m,            !!Placeholder
    vol::Float64 = 2;   # m^3,          !!Placeholder
    m::Float64 = 0.1    # kg            !!Placeholder
end
mutable struct balloonState
    pos::Vector{Float64} = [0.0,0.0,0.0]    # inertial position, m
    vel::Vector{Float64} = [0.0,0.0,0.0]    # inertial velocity, m/s
end
mutable struct balloonStateDeriv
    vel::Vector{Float64}                    # inertial velocity, m/s
    acc::Vector{Float64}                    # inertial acceleration m/s^2
end
function balloonstatederiv(x::balloonState, prm::balloonParams)
    # buoyancy force 
    Fb = (prm.ρa-prm.ρh)*g*prm.vol                  # 
    # drag force
    Vrel = x.vel - windField(x.pos); 
    Fd = prm.Cd * 0.5 * prm.ρa * Vrel.^2 * prm.R;   # F = 1/2 * ρa * Vrel^2  * S * Cd
    a = (Fb+F)/m;                                   # a = ∑F/m
    xdot = balloonStateDeriv(x.vel,a);
    return xdot
end
begin
    println("Simulating...");
    history = Vector{balloonState}()
    state = balloonState();
    for i in 1:100
        #  integrate
        state = RK4_int(0.1, state,balloonStateDeriv);
        push!(history, state)
    end
end