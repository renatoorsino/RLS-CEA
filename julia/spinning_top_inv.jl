using LinearAlgebra
using Random
using DifferentialEquations
using Plots


## -- SYSTEM PARAMETERS -- ##

n = 6        # no. of coordinates
c = 3        # no. of constraint equations


## -- PHYSICAL PARAMETERS -- ##

I = 1.0                                                 # [kg*m^2] - moment of inertia (transverse axes)
J = 2 * I                                               # [kg*m^2] - moment of inertia (inertia symmetry axis)
g = 10.0                                                # [m/s^2] - acceleration of gravity
mgd = 0.2                                               # [N*m] - weight * distance to the center of mass
d = 1.0                                                 # [m] - distance to the center of mass
m = mgd / (g * d)                                       # [kg] - mass

q0 = pi / 6                                             # [rad] - initial nutation angle
ws = 1                                                  # [rad/s] - spinning rate 
D = (J * ws)^2 + 4 * mgd * (J - I) * cos(q0)            # [kg^2*m^4/s^2] - delta
wps = (-(J * ws) - sqrt(D)) / (2 * (J - I) * cos(q0))   # [rad/s] - stationary precession rate
wp = 0.00 * wps                                         # [rad/s] - initial precession rate

K3 = (J * (ws + wp * cos(q0)))                          # [kg*m^2/s] - angular momentum (inertia symmetry axis)
Kz = (I * wp * (sin(q0))^2 + K3 * cos(q0))              # [kg*m^2/s] - angular momentum (vertical axis)

bg = [2 * 1 * 10.0, 10.0^2]                             # Baumgarte constants (ζ = 1, ω = 10)

## -- EQUATION OF MOTION -- ##    

dy = zeros(2 * n)
ac = zeros(n)

M = diagm([m, m, m, I - m * d^2, I - m * d^2, J])
f = [0, 0, -m * g, 0, 0, 0]

h0 = zeros(c)
h1 = zeros(c)
A = diagm(c, n, ones(c))
b = zeros(c)


function constraint_update!(t, q, v)
    # update A 
    A[1, 4] = -d * sin(q[4])
    A[1, 5] = -d * cos(q[5]) * cos(q[4])
    A[2, 4] = d * cos(q[4])
    A[2, 5] = -d * cos(q[5]) * sin(q[4])
    A[3, 5] = d * sin(q[5])

    # update b
    b[1] = -d * (-v[4] * v[5] * sin(q[4]) / tan(q[5]) + v[4]^2 * cos(q[4]) / sin(q[5]) + v[5]^2 * sin(q[5]) * cos(q[4]))
    b[2] = -d * (v[4] * v[5] * cos(q[4]) / tan(q[5]) + v[4]^2 * sin(q[4]) / sin(q[5]) + v[5]^2 * sin(q[5]) * sin(q[4]))
    b[3] = -d * v[5]^2 * cos(q[5])

    # update h0
    h0[1] = q[1] - d * cos(q[4]) * sin(q[5])
    h0[2] = q[2] - d * sin(q[4]) * sin(q[5])
    h0[3] = q[3] - d * cos(q[5])

    # update h1
    h1[1] = v[1] + A[1, 4] * v[4] + A[1, 5] * v[5]
    h1[2] = v[2] + A[2, 4] * v[4] + A[2, 5] * v[5]
    h1[3] = v[3] + A[3, 5] * v[5]
end


function spinning_top!(dy, y, p, t)
    I, J, d = p

    dy[1:3] = copy(y[7:9])
    dy[4] = -y[10] / sin(y[5])
    dy[5] = y[11]
    dy[6] = y[12] - dy[4] * cos(y[5])

    # update f 
    Ho = (I * dy[4] * cos(y[5]) - J * y[12])
    f[4] = Ho * y[11]
    f[5] = -Ho * y[10]

    # update constraint equations
    constraint_update!(t, y[1:n], y[(n+1):(2*n)])

    ac[1:6] = copy(copy((hcat([M; A], [A'; zeros(c, c)])\[f; b])[1:6]))
    dy[7:12] = copy(copy((hcat([M; A], [A'; zeros(c, c)])\[f; b - bg[1] * h1 - bg[2] * h0])[1:6]))
end


## -- INTEGRATOR SETUP - TIME VECTOR & INITIAL CONDITIONS -- ##

t0 = 0.0                        # [s] - initial time
t1 = 1.0e4 # 250.0              # [s] - final time
# dt = 10 # 0.05                  # [s] - sampling time

y0 = [d * sin(q0), 0, d * cos(q0), 0, q0, 0, 0, wp * d * sin(q0), 0, -wp * sin(q0), 0, ws + wp * cos(q0)]           # initial conditions

prob = ODEProblem(spinning_top!, y0, (t0, t1), (I, J, d))

begin
    it = init(prob, DP8(), reltol=1e-30)
    t = [0.0]
    sol = [y0]
    n0 = [0.0]
    n1 = [0.0]
    n2 = [0.0]
    while it.t < t1
        step!(it)
        push!(n2, norm(b - A * ac))
        q = copy(it.u[1:n])
        v = copy(it.u[(n+1):(2*n)])
        constraint_update!(t, q, v)
        push!(sol, copy([q; v]))
        push!(t, it.t)
        push!(n0, norm(h0))
        push!(n1, norm(h1))
        reinit!(it, copy(it.u), t0=it.t)
    end
end


# sol = solve(prob, saveat=dt, reltol=1e-30)
# t = sol.t
# Y = copy(Array(sol))

## -- PLOT -- ##

Y = hcat(sol...)

K3 = @. J * (Y[12, :])
Kz = @. (-I * Y[10, :] * sin(Y[5, :]) + K3 * cos(Y[5, :]))
ME = @. (0.5 * (I * Y[10, :]^2 + I * Y[11, :]^2 + J * Y[12, :]^2) + mgd * cos(Y[5, :]))

p_inv = Dict()

# p_inv["(θ,ω)"] = plot(
#     sol, vars=(5, 11), color="indianred1",
#     xlabel="θ [rad]", ylabel="ω [rad/s]", label=false
# )

p_inv["xyz"] = plot(
    Y[1, :], Y[2, :], Y[3, :], color="steelblue1",
    xlabel="x/d", ylabel="y/d", zlabel="z/d",
    xlims=[-1, 1], ylims=[-1, 1], zlims=[0, 1],
    label=false
)

# begin
#     p_inv["θ"] = plot(
#         sol, vars=(0, 5), color="gray60",
#         xlabel="t [s]", ylabel="θ [rad]", label=false
#     )
#     p_inv["ϕ"] = plot(
#         sol, vars=(0, 4), color="coral3",
#         xlabel="t [s]", ylabel="ϕ [rad]", label=false
#     )
#     p_inv["ψ"] = plot(
#         sol, vars=(0, 6), color="deepskyblue4",
#         xlabel="t [s]", ylabel="ψ [rad]", label=false
#     )
#     p_inv["s"] = plot(
#         p_inv["θ"], p_inv["ϕ"], p_inv["ψ"],
#         layout=(3, 1)
#     )
# end

begin
    p_inv["E"] = plot(
        t, (ME .- ME[1]), color="gray60",
        xlabel="t [s]", ylabel="ΔE [J]", label=false
    )
    p_inv["K3"] = plot(
        t, K3 .- K3[1], color="coral3",
        xlabel="t [s]", ylabel="ΔK_3 [kg m²/s]", label=false
    )
    p_inv["Kz"] = plot(
        t, Kz .- Kz[1], color="deepskyblue4",
        xlabel="t [s]", ylabel="ΔK_z [kg m²/s]", label=false
    )
    p_inv["s"] = plot(
        p_inv["E"], p_inv["K3"], p_inv["Kz"],
        layout=(3, 1)
    )
end


begin
    p_inv["n0"] = plot(
        t, n0, color="gray60",
        xlabel="t [s]", ylabel="n0", label=false
    )
    p_inv["n1"] = plot(
        t, n1, color="coral3",
        xlabel="t [s]", ylabel="n1", label=false
    )
    p_inv["n2"] = plot(
        t, n2, color="coral3",
        xlabel="t [s]", ylabel="n2", label=false
    )
    p_inv["nn"] = plot(
        p_inv["n0"], p_inv["n1"], p_inv["n2"],
        layout=(3, 1)
    )
end