using DynamicPolynomials, MosekTools, MomentOpt, DifferentialEquations, PyPlot, LaTeXStrings, MarkovBounds
DE = DifferentialEquations
@polyvar(x[1:2])
@polyvar(u[1:1])
@polyvar(t)

μ = 1.0
f = [1/μ*x[2]; μ*(x[2] - 1/3*x[2]^3 - x[1]) + x[1]*u[1]]
g = [0; 0.1*x[1]]
σ = g*g'

Tf = 5.0
x0 = [0.5; 0.5]
U = @set(u[1] <= 1 && u[1] >= -1)

lagrange = x[1]^2 + x[2]^2/10 + u[1]^2/10

## data
order = 6
nt = 10
trange = range(Tf/nt, stop = Tf, length=nt)

## solve
MP = DiffusionProcess(x, f, σ, FullSpace())
CP = ControlProcess(MP, Tf, u, t, U, nothing, nothing, LagrangeMayer(lagrange, 0), 0)

bnds, V = optimal_control(CP, x0, order, trange, Mosek.Optimizer; value_function = true)

N = 10
urange = -1.0:0.05:1.0
function control_policy(z,s)
    one_step_MPC = MarkovBounds.extended_inf_generator(MP, V(z,s)[1], t) + lagrange
    u = urange[findmin([one_step_MPC(z..., u, s) for u in urange])[2]]
    return u
end
drift(x,p,t) = [f[1](x[2]); f[2](x[1],x[2],control_policy(x[1:2], t)); lagrange(x[1],x[2],control_policy(x[1:2], t))]
drift_uncontrolled(x,p,t) = [f[1](x[2]); f[2](x[1],x[2],0); lagrange(x[1],x[2],0)]
diffusion(x,p,t) = [0; g[2](x[1]); 0]
sde = DE.SDEProblem(drift, diffusion, [x0..., 0], (0.0,Tf))
sde_uncontrolled = DE.SDEProblem(drift_uncontrolled, diffusion, [x0..., 0], (0.0,Tf))
sde_ensemble = DE.EnsembleProblem(sde)
sde_ensemble_uncontrolled = DE.EnsembleProblem(sde_uncontrolled)

trajectories = DE.solve(sde_ensemble, trajectories=N)
trajectories_uncontrolled = DE.solve(sde_ensemble_uncontrolled, trajectories=N)

fig, ax = subplots()
for sol in trajectories_uncontrolled
    ax.plot([u[1] for u in sol.u], [u[2] for u in sol.u], color="blue", linewidth=0.5)
end
for sol in trajectories
    ax.plot([u[1] for u in sol.u], [u[2] for u in sol.u], color="red", linewidth=0.5)
end
ax.set_xlabel(L"x_1")
ax.set_ylabel(L"x_2")

empiric_control_cost = sum(sol.u[end][3] for sol in trajectories)/N
display(fig)
println("Empiric control cost = ", empiric_control_cost)
println("Lower bound on control cost = ", bnds.bounds)
println("Degree of suboptimality = ", (empiric_control_cost-bnds.bounds)/bnds.bounds)
