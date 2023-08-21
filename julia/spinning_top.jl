using LinearAlgebra
using DifferentialEquations
using Plots


## -- PHYSICAL PARAMETERS -- ##

I = 1                                                   # [kg*m^2] - moment of inertia (transverse axes)
J = 2 * I                                               # [kg*m^2] - moment of inertia (inertia symmetry axis)
mgd = 0.2                                               # [N*m] - weight * distance to the center of mass

q0 = pi / 6                                             # [rad] - initial nutation angle
ws = 1                                                  # [rad/s] - spinning rate 
D = (J * ws)^2 + 4 * mgd * (J - I) * cos(q0)            # [kg^2*m^4/s^2] - delta
wps = (-(J * ws) - sqrt(D)) / (2 * (J - I) * cos(q0))   # [rad/s] - stationary precession rate
wp = 0.00 * wps                                         # [rad/s] - initial precession rate

K3 = (J * (ws + wp * cos(q0)))                          # [kg*m^2/s] - angular momentum (inertia symmetry axis)
Kz = (I * wp * (sin(q0))^2 + K3 * cos(q0))              # [kg*m^2/s] - angular momentum (vertical axis)


## -- EQUATION OF MOTION -- ##    

function spinning_top!(dy, y, p, t)
    I, J, mgd, Kz, K3 = p
    dy[1] = y[4]
    dy[2] = (Kz - K3 * cos(y[1])) / (I * (sin(y[1]))^2)
    dy[3] = K3 / J - dy[2] * cos(y[1])
    dy[4] = (mgd / I) * sin(y[1]) - (Kz - K3 * cos(y[1])) * (K3 - Kz * cos(y[1])) / (I^2 * (sin(y[1]))^3)
end


## -- INTEGRATOR SETUP - TIME VECTOR & INITIAL CONDITIONS -- ##

t0 = 0.0                        # [s] - initial time
t1 = 1.0e5 # 250.0              # [s] - final time
dt = 10 # 0.05                  # [s] - sampling time

y0 = [q0, 0, 0, 0]              # initial conditions

prob = ODEProblem(spinning_top!, y0, (t0, t1), (I, J, mgd, Kz, K3))
sol = solve(prob, saveat=dt, reltol=1e-30)


## -- PLOT -- ##

t = sol.t
Y = copy(Array(sol))

x = @. sin(Y[1, :]) * cos(Y[2, :])
y = @. sin(Y[1, :]) * sin(Y[2, :])
z = @. cos(Y[1, :])

dϕ = @. (Kz - K3 * cos(Y[1,:])) / (I * (sin(Y[1,:]))^2)
dψ = @. K3 / J - dϕ * cos(Y[1,:])

E = @. (0.5 * (I * (dϕ * sin(Y[1,:]))^2 + I * Y[4,:]^2 + J * (dψ + dϕ * cos(Y[1,:]))^2) + mgd *  cos(Y[1, :]))

pl = Dict()

pl["xyz"] = plot(
    x, y, z, color="steelblue1",
    xlabel="x/d", ylabel="y/d", zlabel="z/d",
    xlims=[-1, 1], ylims=[-1, 1], zlims=[0, 1],
    label=false
)

pl["(θ,ω)"] = plot(
    sol, vars=(1, 4), color="indianred1",
    xlabel="θ [rad]", ylabel="ω [rad/s]", label=false
)       

begin
    pl["θ"] = plot(
        sol, vars=(0, 1), color="gray60",
        xlabel="t [s]", ylabel="θ [rad]", label=false
    )
    pl["ϕ"] = plot(
        sol, vars=(0, 2), color="coral3",
        xlabel="t [s]", ylabel="ϕ [rad]", label=false
    )
    pl["ψ"] = plot(
        sol, vars=(0, 3), color="deepskyblue4",
        xlabel="t [s]", ylabel="ψ [rad]", label=false
    )
    pl["s"] = plot(
        pl["θ"], pl["ϕ"], pl["ψ"], 
        layout = (3,1)
    )    
end

pl["E"] = plot(
    t, (E .- E[1]), color="gray60",
    xlabel="t [s]", ylabel="ΔE [J]", label=false
)