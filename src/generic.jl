using LinearAlgebra, FFMPEG, BenchmarkTools, PyPlot, Optimization, Plots, PyCall
using Base.MathConstants: π
# Runge-Kutta 4th-order integrator
function RK4_int(dt, state, fdot::Function)
    k1 = dt * fdot(state)
    k2 = dt * fdot(state + 0.5 * k1)
    k3 = dt * fdot(state + 0.5 * k2)
    k4 = dt * fdot(state + k3)
    return state + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
end
# define downburst wind field as described in xx et. al.
function WindField(pos::Vector{Float64})
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
# utilities from Himanshu Gupta (thanks~!)
#=
Author: Himanshu Gupta
=#

#Function to get rotation matrix that rotates around the z axis
function GetRotationMatrixFromInertialToV1(yaw)
    rotation_matrix = [
                        cos(yaw) sin(yaw) 0;
                        -sin(yaw) cos(yaw) 0;
                        0 0 1
                        ]
    return rotation_matrix
end

#Function to get rotation matrix that rotates around the y axis
function GetRotationMatrixFromV1ToV2(pitch)
    rotation_matrix = [
                        cos(pitch) 0 -sin(pitch);
                        0 1 0;
                        sin(pitch) 0 cos(pitch)
                        ]
    return rotation_matrix
end

#Function to get rotation matrix that rotates around the x axis
function GetRotationMatrixFromV2ToBody(roll)
    rotation_matrix = [
                        1 0 0;
                        0 cos(roll) sin(roll);
                        0 -sin(roll) cos(roll)
                        ]
    return rotation_matrix
end

#Function to get rotation matrix that rotates around z,y and then x axis
function RotationMatrix321(euler_angles::EulerAngles)

    roll = euler_angles.ϕ
    pitch = euler_angles.θ
    yaw = euler_angles.ψ
    RMFromIntertialToV1 = GetRotationMatrixFromInertialToV1(yaw)
    RMFromV1ToV2 = GetRotationMatrixFromV1ToV2(pitch)
    RMFromV2ToBody = GetRotationMatrixFromV2ToBody(roll)
    #GetRotationMatrixToTransformFromInertialToBody
    RM = RMFromV2ToBody*RMFromV1ToV2*RMFromIntertialToV1
    return RM
end

#Function to transform vector from body frame to inertial frame
function TransformFromBodyToInertial(VectorBody, euler_angles::EulerAngles)
    rotation_matrix = RotationMatrix321(euler_angles)
    v = transpose(rotation_matrix)*VectorBody
    return v
end

#Function to transform vector from inertial frame to body frame
function TransformFromInertialToBody(VectorBody, euler_angles::EulerAngles)
    rotation_matrix = RotationMatrix321(euler_angles)
    v = rotation_matrix*VectorBody
    return v
end

#Function to get Air Relative Vehcile Velocity Vector in Body Frame from Wind Angles
function WindAnglesToAirRelativeVelocityVector(wind_angles::WindAngles)
    Va = wind_angles.Va
    beta = wind_angles.β
    alpha = wind_angles.α
    u = Va*cos(alpha)*cos(beta)
    v = Va*sin(beta)
    w = Va*sin(alpha)*cos(beta)
    AirSpeedInBodyFrame = [u,v,w]  #This is the VelocityBody vector
    return AirSpeedInBodyFrame
end

#Function to get Wind Angles from Air Relative Vehcile Velocity Vector in Body Frame
function AirRelativeVelocityVectorToWindAngles(VelocityBody)
    u = VelocityBody[1]
    v = VelocityBody[2]
    w = VelocityBody[3]
    Va = sqrt( (u*u) + (v*v) + (w*w) )
    beta = asin(v/Va)
    alpha = atan(w/u)
    wind_angles = WindAngles(Va, beta, alpha)
    return wind_angles
end

#Function to get the multiplication matrix for rotational kinematics in ODEs
function GetRotationalKinematicsMatrix(euler_angles::EulerAngles)

    roll = euler_angles.ϕ #phi
    pitch = euler_angles.θ #theta
    yaw = euler_angles.ψ #psi
    multiplication_matrix = [
                        1 sin(roll)*tan(pitch) cos(roll)*tan(pitch);
                        0 cos(roll) -sin(roll);
                        0 sin(roll)*sec(pitch) cos(roll)*sec(pitch)
                        ]
    return multiplication_matrix
end


#Function to get all the gamma values for rotational dynamics
function GetGammaValues(aircraft_parameters::AircraftParameters)
    Ix = aircraft_parameters.Ix
    Iy = aircraft_parameters.Iy
    Iz = aircraft_parameters.Iz
    Ixz = aircraft_parameters.Ixz

    Gamma = (Ix*Iz) - (Ixz^2)
    Gamma1 = (Ixz * (Ix-Iy+Iz))/Gamma
    Gamma2 = ( (Ixz^2) + (Iz*(Iz-Iy)) )/Gamma
    Gamma3 = Iz/Gamma
    Gamma4 = Ixz/Gamma
    Gamma5 = (Iz - Ix)/Iy
    Gamma6 = Ixz/Iy
    Gamma7 = ( (Ixz^2) + (Ix*(Ix-Iy)) )/Gamma
    Gamma8 = Ix/Gamma

    GammaArray = [Gamma1,Gamma2,Gamma3,Gamma4,Gamma5,Gamma6,Gamma7,Gamma8,Gamma]
    return GammaArray
end