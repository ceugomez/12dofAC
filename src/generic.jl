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