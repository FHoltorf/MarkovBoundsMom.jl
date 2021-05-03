using MosekTools, Catalyst, Plots, MarkovBounds
# use Catalyst to define the reaction network under consideration
# in this case, we investigate Schlögel's system
rn = @reaction_network begin
    0.01*2, 2X --> X
    1.0, ∅ --> X
end
X = species(rn)[1]
# specify initial molecular counts in system
x0 = Dict(X => 1)
x_scale =  Dict(X => 10.0)

# compute bounds on mean molecular counts of the product counts at steady state
# for different truncation orders 2:5
orders = 2:2:6
sol = Dict()
for m in orders
    sol[m] = stationary_mean(rn, [X], x0, m, Mosek.Optimizer, auto_scaling = false, scales = x_scale)
end

p = Plots.plot(orders, [sol[m].bounds[X][1] for m in orders], label="lower bound", xlabel = "truncation order", ylabel = "⟨X∞⟩")
plot!(p, orders, [sol[m].bounds[X][2] for m in orders], label="upper bound", legend = :bottomright)
display(p)

# same for variance
sol = Dict()
for m in orders
    sol[m] = stationary_variance(rn, [X], x0, m, Mosek.Optimizer, auto_scaling = false, scales = x_scale)
end
p = Plots.plot(xlabel="truncation order", ylabel=string("Var(X∞)"), legend = :right)
plot!(p, orders, [sol[m].bounds[X] for m in orders], color = :red, label="upper bound")
display(p)

# transient moments
nT = 10
traj = []
Tf_range = [0.0, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 5.0, 7.5,
            10.0, 12.5, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0]
order = 4
for Tf in Tf_range
    trange = range(Tf/nT, stop=Tf, length=nT)
    sol = transient_mean(rn, [X], x0, order, trange, Mosek.Optimizer, scales = x_scale)
    push!(traj, sol)
end
p = Plots.plot(xlabel="time", ylabel="⟨X⟩", legend = :right)
plot!(p, Tf_range, [sol.bounds[X][1] for sol in traj], color = :red, label = "lower bound")
plot!(p, Tf_range, [sol.bounds[X][2] for sol in traj], color = :blue, label = "upper bound")
display(p)

# variance
traj = []
for Tf in Tf_range
    trange = range(Tf/nT, stop=Tf, length=nT)
    sol = transient_variance(rn, [X], x0, order, trange, Mosek.Optimizer, scales =  x_scale)
    push!(traj, sol)
end
p = Plots.plot(xlabel="time", ylabel=string("Var(X)"), legend = :right)
plot!(p, Tf_range, [sol.bounds[X][1] for sol in traj], color = :red, label = "lower bound")
display(p)
