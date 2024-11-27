include("generic.jl");
using Parameters
# State derivative of a balloon in a wind field (to be implemented)
@with_kw mutable struct balloonParams
    # known from swenson et al; 
    # payload mass: 91g 
    # 125 L volume, ~ 0.125 m^3
    gravity::Float64 = 9.81
    ρa::Float64 = 1.225 # kg/m^3        !! Modifiable based on state
    ρh::Float64 = 0.166 # kg/m^3, helium
    Cd::Float64 = 0.01  # drag coeff    !!Placeholder
    r::Float64 = 1;     # m,            !!Placeholder
    vol::Float64 = 0.125;   # m^3,         
    m::Float64 = 0.92    # kg       
end
function balloonstatederiv(x::Vector{Float64})
    # buoyancy force 
    fb = (prm.ρa-prm.ρh)*prm.gravity*prm.vol                  
    Fb = Vector{Float64}([0.0,0.0,fb])                                      # only in z dir
    # drag force
    Vrel = x[4:6] - windField(x[1:3]); 
    Fd = prm.Cd * 0.5 * prm.ρa * norm(Vrel.^2) * prm.r;                     # F = 1/2 * ρa * Vrel^2  * S * Cd
    Fd = -normalize(Vrel)*Fd;                                               # opposes relative motion
    # a = ∑F/m
    a = (Fb+Fd)./prm.m;                                
    # return state derivative
    xdot = [x[4:6]; a]; #vcat
    return xdot
end
function getDensity(z)

    
end
function VisualizeFieldXYZ(history::Matrix{Float64})
    plt = figure();
    title("Balloon Ground Track");
    # NE track
        ax1 = plt.add_subplot(131);
        ax1.plot(history[1,:], history[2,:]);
        xlabel("Easting [m]")
        ylabel("Northing [m]")
        # EU track
        ax2 = plt.add_subplot(132);
        ax2.plot(history[2,:], history[3,:]);
        xlabel("Easting [m]")
        ylabel("Up [m]")
        # 3d plot
        ax3 = plt.add_subplot(133);
        ax3.plot(history[1,:],history[3,:]);
        xlabel("Northing[m]");
        ylabel("Up [m]")
    savefig("./figures/Track")
    show()
    return nothing
end
function VisualizeFieldT(history::Matrix{Float64})
    plt = figure();
    title("Balloon Time History");
    t = LinRange(0, 3600/60, 3601)

    # NE track
        ax1 = plt.add_subplot(131);
        ax1.plot(t, history[1,:]);
        xlabel("Time [s]")
        ylabel("Easting [m]")
        # EU track
        ax2 = plt.add_subplot(132);
        ax2.plot(t, history[2,:]);
        xlabel("Time [s]")
        ylabel("Northing [m]")
        # 3d plot
        ax3 = plt.add_subplot(133);
        ax3.plot(t,history[3,:]);
        xlabel("Time [s]");
        ylabel("Up [m]")
    savefig("./figures/TTrack")
    show()
    return nothing
end
begin

    println("Simulating...");
    # initialize parameter structure
    global prm = balloonParams();
    # initialize history vector
    tstep = 0.1;
    tf = 3600
    history = zeros(Float64, 6,tf+1);

    for i in 1:(tf)
        #  integrate
        stateK = history[:,i];
        stateKp1 = RK4_int(tstep, stateK, balloonstatederiv)
        history[:,i+1] = stateKp1; # integrate
    end
    println("done")
    #VisualizeFieldXYZ(history)
    VisualizeFieldT(history)

end