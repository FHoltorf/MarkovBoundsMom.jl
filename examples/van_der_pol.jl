using MosekTools, DifferentialEquations, PyPlot, LaTeXStrings, MarkovBounds
DE = DifferentialEquations
u = []
@polyvar(x[1:2])
@polyvar(t)

μ = 1/2
f = [x[2]; - x[1] + μ*(1 - x[1]^2) * x[2]]
g = [0; 0.5]
σ = g*g'

Tf = 2.0
x0 = [0.5; 0.5]
U = FullSpace()
X = @set(x[1] >= 0)

## data
order = 8
nt = 20
trange = range(Tf/nt, stop = Tf, length=nt)

## solve
MP = DiffusionProcess(x, f, σ, FullSpace())
CP = ControlProcess(MP, Tf, u, t, U, nothing, nothing, ExitProbability(X), 0)

bnds, V = optimal_control(CP, x0, order, trange, Mosek.Optimizer; value_function = true)


N = 100
drift(x,p,t) = [f[1](x[2]); f[2](x[1],x[2])]
diffusion(x,p,t) = [0; g[2]]
sde = DE.SDEProblem(drift, diffusion, x0, (0.0,Tf))
sde_ensemble = DE.EnsembleProblem(sde)

trajectories = DE.solve(sde_ensemble, trajectories=N, EM(), dt = 0.01)

fig, ax = subplots(figsize=(5,5))
stayed = 0
for sol in trajectories
    max_idx = findfirst(m -> m[1] < 0, sol.u)
    if max_idx == nothing
        global stayed += 1
    end
    max_idx = (max_idx == nothing ? length(sol.u) : max_idx)
    global ln = ax.plot([sol.u[i][1] for i in 1:max_idx], [sol.u[i][2] for i in 1:max_idx], color="blue", linewidth=0.5)
end
ax.scatter([x0[1]], [x0[2]], marker="o", color="black", s=5)
ax.set_xlabel(L"x_1")
ax.set_ylabel(L"x_2")
display(fig)

println("Upper Bound:  P(exit) ≦ ", bnds.bounds)
println("Sampled: P(exit) = ", stayed/N)
