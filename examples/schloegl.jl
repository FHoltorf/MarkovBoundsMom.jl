using MarkovBounds, MosekTools, Catalyst, Plots
# use Catalyst to define the reaction network under consideration
# in this case, we investigate Schlögel's system
rn = @reaction_network begin
    0.15, 2X --> 3X
    1.5e-3, 3X --> 2X
    20.0, ∅ --> X
    2.0, X --> ∅
end
X = species(rn)[1]
# specify initial molecular counts in system
x0 = Dict(X => 20)

# compute bounds on mean molecular counts of the product counts at steady state
# for different truncation orders 2:5
orders = 2:8
sol = mean(rn, [X], x0, orders, Mosek.Optimizer, auto_scaling = false, scales =  Dict(X => 20.0))

pX = plot(orders, [sol[m][1][1] for m in orders], label="lower bound", xlabel = "truncation order", ylabel = "⟨X∞⟩")
plot!(pX, orders, [sol[m][1][2] for m in orders], label="upper bound", legend = :bottomright)
