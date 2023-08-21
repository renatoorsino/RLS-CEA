using LinearAlgebra
using Random
using DifferentialEquations
using DelimitedFiles
using Plots


## -- SYSTEM PARAMETERS -- ##

n = 10        # no. of coordinates
c = 10        # no. of constraint equations


## -- PHYSICAL PARAMETERS -- ##

m = 1.0                       # [kg] - mass of each bar
g = 9.8                       # [m/s^2] - acceleration of gravity
l = 1.0                       # [m] - length of each bar
r = 0.1                       # [m] - square cross section width
k = 1 + (r / l)^2
mo = m * (3 - k) / 12
md = m * (3 + k) / 6
mg = m * g

bg = [2 * 1 * 10.0, 10.0^2]   # Baumgarte constants (ζ = 1, ω = 10)

## -- EQUATION OF MOTION -- ##    

dy = zeros(2 * n)

P = Dict()
z = Dict()
h = Dict()
a = Dict()

s = zeros(c)
e = zeros(c)

M = zeros(n, n)
f = zeros(n)

M[1, 1] = md
M[1, 3] = mo
M[2, 2] = md
M[2, 5] = mo
M[3, 1] = mo
M[3, 3] = md
M[3, 6] = mo
M[4, 4] = md
M[4, 7] = mo
M[5, 2] = mo
M[5, 5] = md
M[5, 8] = mo
M[6, 3] = mo
M[6, 6] = md
M[7, 4] = mo
M[7, 7] = md
M[7, 9] = mo
M[8, 5] = mo
M[8, 8] = md
M[8, 10] = mo
M[9, 7] = mo
M[9, 9] = md
M[10, 8] = mo
M[10, 10] = md

f[4] = -mg
f[7] = -mg
f[9] = -mg

P[0] = inv(M)
a[0] = P[0] * f

h0 = zeros(c)
h1 = zeros(c)
A = zeros(c, n)
b = zeros(c)
i = Vector{Int64}(undef, c)


function constraint_update!(t, y)
    # update A 
    A[1, 1] = y[1]
    A[1, 2] = y[2]
    A[2, 1] = y[1] - y[3]
    A[2, 2] = (-l + y[2] - y[5]) / 2.0 + (l + y[2] - y[5]) / 2.0
    A[2, 3] = -y[1] + y[3]
    A[2, 4] = y[4]
    A[2, 5] = (-l - y[2] + y[5]) / 2.0 + (l - y[2] + y[5]) / 2.0
    A[3, 3] = y[3] - y[6]
    A[3, 4] = y[4] - y[7]
    A[3, 5] = (-l + y[5] - y[8]) / 2.0 + (l + y[5] - y[8]) / 2.0
    A[3, 6] = -y[3] + y[6]
    A[3, 7] = -y[4] + y[7]
    A[3, 8] = (-l - y[5] + y[8]) / 2.0 + (l - y[5] + y[8]) / 2.0
    A[4, 6] = y[6]
    A[4, 7] = y[7] - y[9]
    A[4, 8] = (-l + y[8] - y[10]) / 2.0 + (l + y[8] - y[10]) / 2.0
    A[4, 9] = -y[7] + y[9]
    A[4, 10] = (-l - y[8] + y[10]) / 2.0 + (l - y[8] + y[10]) / 2.0
    A[5, 9] = y[9]
    A[5, 10] = y[10] / 2.0 + (-2 * l + y[10]) / 2.0
    A[6, 1] = y[5] / 2.0
    A[6, 2] = -0.5 * y[3]
    A[6, 3] = -0.5 * y[2]
    A[6, 5] = y[1] / 2.0
    A[7, 7] = (-l + y[10]) / 2.0
    A[7, 8] = -0.5 * y[9]
    A[7, 9] = (l - y[8]) / 2.0
    A[7, 10] = y[7] / 2.0
    A[8, 1] = y[1] * (-y[4] + y[7]) + (y[4] * y[6] - y[3] * y[7]) / 2.0
    A[8, 2] = -0.5 * (y[2] * y[4]) + (y[2] * y[7]) / 2.0 + ((y[2] - y[5]) * y[7] + y[4] * (-y[2] + y[8])) / 2.0
    A[8, 3] = -0.5 * (y[1] * y[7])
    A[8, 4] = -0.5 * (y[1] * y[1]) + (y[1] * y[6]) / 2.0 + (y[2] * (-y[2] + y[8])) / 2.0
    A[8, 5] = -0.5 * (y[2] * y[7])
    A[8, 6] = (y[1] * y[4]) / 2.0
    A[8, 7] = (y[1] * y[1]) / 2.0 - (y[1] * y[3]) / 2.0 + (y[2] * (y[2] - y[5])) / 2.0
    A[8, 8] = (y[2] * y[4]) / 2.0
    A[9, 3] = ((y[7] - y[9]) * y[9] + (l - y[10]) * (-y[8] + y[10])) / 2.0
    A[9, 4] = -0.5 * (y[6] * y[9])
    A[9, 5] = (y[6] * (l - y[10])) / 2.0
    A[9, 6] = (-(y[4] * y[9]) + y[9] * y[9] + (l - y[10]) * (y[5] - y[10])) / 2.0
    A[9, 7] = (y[3] * y[9]) / 2.0
    A[9, 8] = -0.5 * (y[3] * (l - y[10]))
    A[9, 9] = -0.5 * (y[4] * y[6]) + (y[3] * (y[7] - y[9])) / 2.0 - (y[3] * y[9]) / 2.0 + y[6] * y[9]
    A[9, 10] = (y[3] * (l - y[10])) / 2.0 - (y[6] * (l - y[10])) / 2.0 - (y[6] * (y[5] - y[10])) / 2.0 - (y[3] * (-y[8] + y[10])) / 2.0
    A[10, 1] = y[1] * y[6] * y[9] + (-(y[3] * y[6] * y[9]) - y[4] * ((y[7] - y[9]) * y[9] + (l - y[10]) * (-y[8] + y[10]))) / 2.0
    A[10, 2] = (y[2] * y[6] * y[9]) / 2.0 + (y[6] * ((y[2] - y[5]) * y[9] + y[4] * (-l + y[10]))) / 2.0
    A[10, 3] = -0.5 * (y[1] * y[6] * y[9])
    A[10, 4] = (y[2] * y[6] * (-l + y[10])) / 2.0 - (y[1] * ((y[7] - y[9]) * y[9] + (l - y[10]) * (-y[8] + y[10]))) / 2.0
    A[10, 5] = -0.5 * (y[2] * y[6] * y[9])
    A[10, 6] = (y[1] * y[1] * y[9]) / 2.0 - (y[1] * y[3] * y[9]) / 2.0 + (y[2] * ((y[2] - y[5]) * y[9] + y[4] * (-l + y[10]))) / 2.0
    A[10, 7] = -0.5 * (y[1] * y[4] * y[9])
    A[10, 8] = (y[1] * y[4] * (l - y[10])) / 2.0
    A[10, 9] = (y[1] * y[1] * y[6]) / 2.0 - (y[1] * y[3] * y[6]) / 2.0 + (y[2] * (y[2] - y[5]) * y[6]) / 2.0 - (y[1] * y[4] * (y[7] - y[9])) / 2.0 + (y[1] * y[4] * y[9]) / 2.0
    A[10, 10] = (y[2] * y[4] * y[6]) / 2.0 - (y[1] * y[4] * (l - y[10])) / 2.0 + (y[1] * y[4] * (-y[8] + y[10])) / 2.0

    # update b
    b[1] = -(y[11] * y[11]) - y[12] * y[12]
    b[2] = -((y[11] - y[13]) * (y[11] - y[13])) - y[14] * y[14] - (y[12] - y[15]) * (y[12] - y[15])
    b[3] = -((y[13] - y[16]) * (y[13] - y[16])) - (y[14] - y[17]) * (y[14] - y[17]) - (y[15] - y[18]) * (y[15] - y[18])
    b[4] = -(y[16] * y[16]) - (y[17] - y[19]) * (y[17] - y[19]) - (y[18] - y[20]) * (y[18] - y[20])
    b[5] = -(y[19] * y[19]) - y[20] * y[20]
    b[6] = y[12] * y[13] - y[11] * y[15]
    b[7] = y[18] * y[19] - y[17] * y[20]
    b[8] = -((-y[4] + y[7]) * (y[11] * y[11])) - y[1] * y[14] * y[16] + y[1] * y[13] * y[17] - y[2] * (y[12] - y[15]) * y[17] - 2 * y[1] * y[11] * (-y[14] + y[17]) - y[11] * (-(y[7] * y[13]) + y[6] * y[14] + y[4] * y[16] - y[3] * y[17]) - y[2] * y[14] * (-y[12] + y[18]) - y[12] * ((-y[2] + y[8]) * y[14] + y[7] * (y[12] - y[15]) + (y[2] - y[5]) * y[17] + y[4] * (-y[12] + y[18]))
    b[9] = y[6] * y[14] * y[19] - y[3] * (y[17] - y[19]) * y[19] - y[6] * (y[19] * y[19]) - y[3] * (y[18] - y[20]) * y[20] - y[6] * y[20] * (-y[15] + y[20]) - y[16] * ((l - y[10]) * y[15] - y[9] * (y[14] - 2 * y[19]) - y[4] * y[19] - (l + y[5] - 2 * y[10]) * y[20]) - y[13] * ((-l + y[10]) * y[18] + y[9] * (y[17] - 2 * y[19]) + y[7] * y[19] + (l + y[8] - 2 * y[10]) * y[20])
    b[10] = -(y[6] * y[9] * (y[11] * y[11])) - 2 * y[1] * y[9] * y[11] * y[16] - ((y[2] - y[5]) * y[9] + y[4] * (-l + y[10])) * y[12] * y[16] + y[1] * y[9] * y[13] * y[16] - 2 * y[1] * y[6] * y[11] * y[19] + y[1] * y[6] * y[13] * y[19] - y[2] * y[6] * (y[12] - y[15]) * y[19] - y[1] * y[1] * y[16] * y[19] + y[1] * y[3] * y[16] * y[19] - y[1] * y[4] * y[19] * (-y[17] + y[19]) - y[2] * y[6] * y[14] * y[20] - y[1] * y[4] * y[20] * (-y[18] + y[20]) - y[6] * y[12] * ((-l + y[10]) * y[14] + y[9] * (y[12] - y[15]) + (y[2] - y[5]) * y[19] + y[4] * y[20]) - y[2] * y[16] * ((-l + y[10]) * y[14] + y[9] * (y[12] - y[15]) + (y[2] - y[5]) * y[19] + y[4] * y[20]) + y[1] * y[14] * ((-l + y[10]) * y[18] + y[9] * (y[17] - 2 * y[19]) + y[7] * y[19] + (l + y[8] - 2 * y[10]) * y[20]) + y[11] * (y[6] * y[9] * y[13] + ((y[7] - y[9]) * y[9] + (l - y[10]) * (-y[8] + y[10])) * y[14] + y[3] * y[9] * y[16] + y[3] * y[6] * y[19] + y[4] * ((-l + y[10]) * y[18] + y[9] * (y[17] - 2 * y[19]) + y[7] * y[19] + (l + y[8] - 2 * y[10]) * y[20]))

    # update h0
    h0[1] = (-(l * l) + y[1] * y[1] + y[2] * y[2]) / 2.0
    h0[2] = ((y[1] - y[3]) * (y[1] - y[3]) + y[4] * y[4] - (l + y[2] - y[5]) * (l - y[2] + y[5])) / 2.0
    h0[3] = ((y[3] - y[6]) * (y[3] - y[6]) + (y[4] - y[7]) * (y[4] - y[7]) - (l + y[5] - y[8]) * (l - y[5] + y[8])) / 2.0
    h0[4] = (y[6] * y[6] + (y[7] - y[9]) * (y[7] - y[9]) - (l + y[8] - y[10]) * (l - y[8] + y[10])) / 2.0
    h0[5] = (y[9] * y[9] + y[10] * (-2 * l + y[10])) / 2.0
    h0[6] = (-(y[2] * y[3]) + y[1] * y[5]) / 2.0
    h0[7] = ((l - y[8]) * y[9] + y[7] * (-l + y[10])) / 2.0
    h0[8] = (y[1] * y[1] * (-y[4] + y[7]) + y[1] * (y[4] * y[6] - y[3] * y[7]) + y[2] * ((y[2] - y[5]) * y[7] + y[4] * (-y[2] + y[8]))) / 2.0
    h0[9] = (y[6] * (-(y[4] * y[9]) + y[9] * y[9] + (l - y[10]) * (y[5] - y[10])) + y[3] * ((y[7] - y[9]) * y[9] + (l - y[10]) * (-y[8] + y[10]))) / 2.0
    h0[10] = (y[1] * y[1] * y[6] * y[9] + y[2] * y[6] * ((y[2] - y[5]) * y[9] + y[4] * (-l + y[10])) - y[1] * (y[3] * y[6] * y[9] + y[4] * ((y[7] - y[9]) * y[9] + (l - y[10]) * (-y[8] + y[10])))) / 2.0

    # update h1
    h1[1] = y[1] * y[11] + y[2] * y[12]
    h1[2] = y[2] * y[12] - y[5] * y[12] + y[1] * (y[11] - y[13]) + y[3] * (-y[11] + y[13]) + y[4] * y[14] - y[2] * y[15] + y[5] * y[15]
    h1[3] = (y[3] - y[6]) * (y[13] - y[16]) + (y[4] - y[7]) * (y[14] - y[17]) + ((l + y[5] - y[8]) * (y[15] - y[18])) / 2.0 - ((l - y[5] + y[8]) * (y[15] - y[18])) / 2.0
    h1[4] = y[6] * y[16] - y[9] * y[17] + y[8] * y[18] - y[10] * y[18] + y[7] * (y[17] - y[19]) + y[9] * y[19] - y[8] * y[20] + y[10] * y[20]
    h1[5] = y[9] * y[19] + (-l + y[10]) * y[20]
    h1[6] = (y[5] * y[11] - y[3] * y[12] - y[2] * y[13] + y[1] * y[15]) / 2.0
    h1[7] = ((-l + y[10]) * y[17] - y[9] * y[18] + (l - y[8]) * y[19] + y[7] * y[20]) / 2.0
    h1[8] = (2 * y[1] * (-y[4] + y[7]) * y[11] + (y[4] * y[6] - y[3] * y[7]) * y[11] + ((y[2] - y[5]) * y[7] + y[4] * (-y[2] + y[8])) * y[12] + y[1] * y[1] * (-y[14] + y[17]) + y[1] * (-(y[7] * y[13]) + y[6] * y[14] + y[4] * y[16] - y[3] * y[17]) + y[2] * ((-y[2] + y[8]) * y[14] + y[7] * (y[12] - y[15]) + (y[2] - y[5]) * y[17] + y[4] * (-y[12] + y[18]))) / 2.0
    h1[9] = (((y[7] - y[9]) * y[9] + (l - y[10]) * (-y[8] + y[10])) * y[13] + (-(y[4] * y[9]) + y[9] * y[9] + (l - y[10]) * (y[5] - y[10])) * y[16] + y[6] * (-(y[9] * y[14]) - y[4] * y[19] + 2 * y[9] * y[19] + (l - y[10]) * (y[15] - y[20]) - (y[5] - y[10]) * y[20]) + y[3] * (-(l * y[18]) + y[10] * y[18] + y[9] * (y[17] - 2 * y[19]) + y[7] * y[19] + l * y[20] + y[8] * y[20] - 2 * y[10] * y[20])) / 2.0
    h1[10] = (2 * y[1] * y[6] * y[9] * y[11] - (y[3] * y[6] * y[9] + y[4] * ((y[7] - y[9]) * y[9] + (l - y[10]) * (-y[8] + y[10]))) * y[11] + y[6] * ((y[2] - y[5]) * y[9] + y[4] * (-l + y[10])) * y[12] + y[1] * y[1] * y[9] * y[16] + y[2] * ((y[2] - y[5]) * y[9] + y[4] * (-l + y[10])) * y[16] + y[1] * y[1] * y[6] * y[19] + y[2] * y[6] * ((-l + y[10]) * y[14] + y[9] * (y[12] - y[15]) + (y[2] - y[5]) * y[19] + y[4] * y[20]) - y[1] * (y[6] * y[9] * y[13] + ((y[7] - y[9]) * y[9] + (l - y[10]) * (-y[8] + y[10])) * y[14] + y[3] * y[9] * y[16] + y[3] * y[6] * y[19] + y[4] * (-(l * y[18]) + y[10] * y[18] + y[9] * (y[17] - 2 * y[19]) + y[7] * y[19] + l * y[20] + y[8] * y[20] - 2 * y[10] * y[20]))) / 2.0

end


mechanical_energy(y) = (2 * mg * (y[4] + y[7] + y[9]) + md * (y[11] * y[11]) + md * (y[12] * y[12]) + 2 * mo * y[11] * y[13] + md * (y[13] * y[13]) + md * (y[14] * y[14]) + 2 * mo * y[12] * y[15] + md * (y[15] * y[15]) + 2 * mo * y[13] * y[16] + md * (y[16] * y[16]) + 2 * mo * y[14] * y[17] + md * (y[17] * y[17]) + 2 * mo * y[15] * y[18] + md * (y[18] * y[18]) + 2 * mo * y[17] * y[19] + md * (y[19] * y[19]) + 2 * mo * y[18] * y[20] + md * (y[20] * y[20])) / 2.0


function bricard!(dy, y, p, t)
    l = p

    dy[1:n] = copy(y[(n+1):(2*n)])

    # update constraint equations
    constraint_update!(t, y)

    # shuffle constraint equations
    i[1:c] = copy(shuffle(1:c))

    # recursive constraint enforcement algorithm
    for r in 0:(c-1)
        h[r+1] = A[i[r+1], :]
        z[r+1] = P[r] * h[r+1]
        s[r+1] = h[r+1]' * z[r+1]
        e[r+1] = (b[i[r+1]] - h[r+1]' * a[r] - bg[1] * h1[i[r+1]] - bg[2] * h0[i[r+1]])
        if abs(s[r+1]) > 8e-17
            P[r+1] = P[r] - z[r+1] * (z[r+1]' / s[r+1])
            a[r+1] = a[r] + z[r+1] * (e[r+1] / s[r+1])
        else
            P[r+1] = P[r]
            a[r+1] = a[r]
        end
    end

    dy[(n+1):(2*n)] = copy(a[c])
end


## -- INTEGRATOR SETUP - TIME VECTOR & INITIAL CONDITIONS -- ##

t0 = 0.0                        # [s] - initial time
t1 = 10.0                       # [s] - final time

# initial conditions
y0 = zeros(2 * n)
y0[1:n] = [l, 0, l, -l, 0, l, -l, l, -l, l]

prob = ODEProblem(bricard!, y0, (t0, t1), (l))

begin
    it = init(prob, DP8(), reltol=1e-30, abstol=1e-16)
    t = [0.0]
    sol = [y0]
    n0 = [0.0]
    n1 = [0.0]
    n2 = [0.0]
    ME = [mechanical_energy(y0)]
    while it.t < t1
        step!(it)
        q = copy(it.u[1:n])
        v = copy(it.u[(n+1):(2*n)])
        constraint_update!(t, [q; v])
        for r in 0:(c-1)
            h[r+1] = A[i[r+1], :]
            z[r+1] = P[r] * h[r+1]
            s[r+1] = h[r+1]' * z[r+1]
            e[r+1] = (b[i[r+1]] - h[r+1]' * a[r])
            if abs(s[r+1]) > 8e-17
                P[r+1] = P[r] - z[r+1] * (z[r+1]' / s[r+1])
                v = copy(v - z[r+1] * (h1[i[r+1]] / s[r+1]))
                q = copy(q - z[r+1] * (h0[i[r+1]] / s[r+1]))
                a[r+1] = a[r] + z[r+1] * (e[r+1] / s[r+1])
            else
                P[r+1] = P[r]
                a[r+1] = a[r]
            end
        end
        constraint_update!(t, [q; v])
        push!(sol, copy([q; v]))
        if Bool(floor(it.t) - floor(t[end]))
            print("*")
        end
        push!(t, it.t)
        push!(n0, norm(h0))
        push!(n1, norm(h1))
        push!(n2, norm(b - A * a[c]))
        push!(ME, mechanical_energy([q; v]))
        reinit!(it, copy([q; v]), t0=it.t)
    end
end

# sol = solve(prob, saveat=dt, reltol=1e-30)
# t = sol.t
# Y = copy(Array(sol))

## -- PLOT -- ##

Y = hcat(sol...)
writedlm("bricard_ce.txt", hcat(t, Y[3, :], Y[4, :], Y[5, :], (ME .- ME[1]), n0, n1, n2), "\t")

p_ce = Dict()

p_ce["xyz"] = plot(
    Y[5, :], Y[3, :], Y[4, :], color="steelblue1",
    xlabel="z_2", ylabel="x_2", zlabel="y_2",
    # xlims=[-1, 1], ylims=[-1, 1], zlims=[0, 1],
    label=false
)

begin
    p_ce["x"] = plot(
        t, Y[3, :], color="gray60",
        xlabel="t [s]", ylabel="x_2 [m]", label=false
    )
    p_ce["y"] = plot(
        t, Y[4, :], color="coral3",
        xlabel="t [s]", ylabel="y_2 [m]", label=false
    )
    p_ce["z"] = plot(
        t, Y[5, :], color="deepskyblue4",
        xlabel="t [s]", ylabel="z_2 [m]", label=false
    )
    p_ce["s"] = plot(
        p_ce["x"], p_ce["y"], p_ce["z"],
        layout=(3, 1)
    )
end

p_ce["E"] = plot(
    t, (ME .- ME[1]), color="gray60", 
    xlabel="t [s]", ylabel="ΔE [J]", label=false
)

begin
    p_ce["n0"] = plot(
        t, n0, color="gray60",
        xlabel="t [s]", ylabel="n0", label=false
    )
    p_ce["n1"] = plot(
        t, n1, color="coral3",
        xlabel="t [s]", ylabel="n1", label=false
    )
    p_ce["n2"] = plot(
        t, n2, color="coral3",
        xlabel="t [s]", ylabel="n2", label=false
    )
    p_ce["nn"] = plot(
        p_ce["n0"], p_ce["n1"], p_ce["n2"],
        layout=(3, 1)
    )
end