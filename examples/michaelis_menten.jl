using Catalyst, LaTeXStrings, MosekTools, MarkovBounds
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
ps = []
for s in [S, P]
    p = Plots.plot(xlabel="truncation order", ylabel=string("⟨",s.f.name,"∞⟩"), legend = :right)
    plot!(p, orders, [sol[m].bounds[s][1] for m in orders], color = :red, label="lower bound")
    plot!(p, orders, [sol[m].bounds[s][2] for m in orders], color = :blue, label="upper bound")
    push!(ps, p)
end
Plots.plot(ps...)

# same for variance
sol = Dict()
for m in orders
    sol[m] = stationary_variance(rn, [S,P], x0, m, Mosek.Optimizer, auto_scaling = false, scales = x_scales)
end
ps = []
for s in [S, P]
    p = Plots.plot(xlabel="truncation order", ylabel=string("Var(",s.f.name,"∞)"), legend = :right)
    plot!(p, orders, [sol[m].bounds[s] for m in orders], color = :blue, label="upper bound")
    push!(ps, p)
end
Plots.plot(ps...)

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

ps = []
for s in [S,P]
    plt = Plots.plot(xlabel="time", ylabel=string("⟨", s.f.name,"⟩"), legend = :right)
    plot!(plt, Tf_range, [sol.bounds[s][1] for sol in traj], color = :red, label = "lower bound")
    plot!(plt, Tf_range, [sol.bounds[s][2] for sol in traj], color = :blue, label = "upper bound")
    push!(ps, plt)
end
Plots.plot(ps...)

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
    plt = Plots.plot(xlabel="time", ylabel=string("Var(",s,")"), legend = :right)
    plot!(plt, Tf_range, [sol.bounds[s][1] for sol in traj], color = :red, label = "lower bound")
    push!(ps, plt)
end
Plots.plot(ps...)


## volume of confidence ellipsoid
orders = 2:8
sol = Dict()
for m in orders
    sol[m] = MarkovBounds.stationary_confidence_ellipsoid(rn, x0, [S,P], m, Mosek.Optimizer; scales = x_scales)
end

plt = plot(orders, [sol[m].bounds[[S,P]] for m in orders], xlabel="order", ylabel="Volume of S-P confidence ellipsoids")
