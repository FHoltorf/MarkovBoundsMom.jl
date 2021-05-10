using LaTeXStrings, PyPlot, MosekTools, MarkovBounds
# use Catalyst to define the reaction network under consideration
# in this case, we investigate simple Michaelis-Menten kinetics
rn = @reaction_network begin
    10.0, S + E --> SE
    0.2, SE --> S + E
    4.0, SE --> P + E
    1.3, P --> S
end
S, E, SE, P = species(rn)
# specify initial molecular counts in system
x0 = Dict(S => 10, E => 10, SE => 0, P => 0)
x_scales = Dict(S => 5, E => 5, SE => 5, P => 5)
# compute bounds on mean molecular counts of the product counts at steady state
# for different truncation orders
orders = 2:6
sol = Dict()
for m in orders
    sol[m] = stationary_mean(rn, [S,P], x0, m, Mosek.Optimizer, auto_scaling = false, scales = x_scales)
end


for s in [S, P]
    fig, ax = subplots()
    ax.plot(orders, [sol[m].bounds[s][1] for m in orders], marker = "*", color="blue", label="lower bound")
    ax.plot(orders, [sol[m].bounds[s][2] for m in orders], marker = "*", color="red", label="upper bound")
    ax.set_xlabel("truncation order")
    ax.set_ylabel(string(L"⟨",s.f.name,L"_{\infty}⟩"))
    ax.legend()
    display(fig)
end



# same for variance
sol = Dict()
for m in orders
    sol[m] = stationary_variance(rn, [S,P], x0, m, Mosek.Optimizer, auto_scaling = false, scales = x_scales)
end

for s in [S, P]
    fig, ax = subplots()
    ax.plot(orders, [sol[m].bounds[s] for m in orders], marker = "*", color="red", label="upper bound")
    ax.set_xlabel("truncation order")
    ax.set_ylabel(string(L"Var(",s.f.name,L"_{\infty})"))
    ax.legend()
    display(fig)
end

# transient moments
nT = 5
traj = []
Tf_range = [0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3,
            0.4, 0.5, 0.625, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 5.0]
order = 4
for Tf in Tf_range
    trange = range(Tf/nT, stop=Tf, length=nT)
    sol = transient_mean(rn, [S,P], x0, order, trange, Mosek.Optimizer, scales = x_scales)
    push!(traj, sol)
end

for s in [S,P]
    fig, ax = subplots()
    ax.plot(Tf_range, [sol.bounds[s][1] for sol in traj], color="blue", label = "lower bound")
    ax.plot(Tf_range, [sol.bounds[s][2] for sol in traj], color="red", label = "upper bound")
    ax.set_xlabel("time")
    ax.set_ylabel(string(L"⟨", s.f.name,L"⟩"))
    ax.legend()
    display(fig)
end

# variance
nT = 5
traj = []
Tf_range = [0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3,
            0.4, 0.5, 0.625, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 5.0]
order = 4
for Tf in Tf_range
    trange = range(Tf/nT, stop=Tf, length=nT)
    sol = transient_variance(rn, [P], x0, order, trange, Mosek.Optimizer, scales =  x_scales)
    push!(traj, sol)
end

ps = []
for s in [P]
    fig, ax = subplots()
    ax.plot(Tf_range, [sol.bounds[s][1] for sol in traj], color="red", label = "upper bound")
    ax.set_xlabel("time")
    ax.set_ylabel(string(L"Var(", s.f.name,L")"))
    ax.legend()
    display(fig)
end

## volume of confidence ellipsoid
orders = 2:8
sol = Dict()
for m in orders
    sol[m] = stationary_confidence_ellipsoid(rn, x0, [S,P], m, Mosek.Optimizer; scales = x_scales)
end

fig, ax = subplots()
ax.plot(orders, [sol[m].bounds[[S,P]] for m in orders], marker="*", color="red", label = "upper bound")
ax.set_xlabel("order")
ax.set_ylabel("Volume of S-P confidence ellipsoids")
ax.legend()
display(fig)
