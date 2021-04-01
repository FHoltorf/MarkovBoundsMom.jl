using MarkovBounds, MosekTools, Catalyst, Plots
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

# compute bounds on mean molecular counts of the product counts at steady state
# for different truncation orders 2:5
orders = 2:6
sol = mean(rn, [S,P], x0, orders, Mosek.Optimizer, auto_scaling = false, scales =  Dict(S => 5, E => 5, SE => 5, P => 5))

pS = plot(orders, [sol[m][1][1] for m in orders], label="lower bound", xlabel = "truncation order", ylabel = "⟨S∞⟩")
plot!(pS, orders, [sol[m][1][2] for m in orders], label="upper bound")
pP = plot(orders, [sol[m][2][1] for m in orders], label="lower bound", xlabel = "truncation order", ylabel = "⟨P∞⟩")
plot!(pP, orders, [sol[m][2][2] for m in orders], label="upper bound", legend = :bottomright)
plot(pS, pP)
