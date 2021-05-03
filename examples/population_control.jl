using DynamicPolynomials, MosekTools, MomentOpt, DifferentialEquations, PyPlot, LaTeXStrings, MarkovBounds
DE = DifferentialEquations
@polyvar(x[1:2])
@polyvar(u[1:1])
@polyvar(t)
γ = [1, 2, 1, 2, 0.25*0.1]

f = [γ[1] * x[1] - γ[2] * x[1] * x[2] ;
     γ[4] * x[1] * x[2] - γ[3] * x[2] - x[2]*u[1]]

g = [γ[5]*x[1]; 0]
σ = polynomial.(g*g')

x0 = [1.0; 0.25]

Tf = 10.0
x_target = [0.75, 0.5]
u_target = 0.5

lagrange = (x[1] - x_target[1])^2 + (x[2] - x_target[2])^2/10 + (u[1] - u_target)^2/10

## data
order = 4
nt = 10
trange = range(Tf/nt, stop = Tf, length=nt)

## solve
MP = DiffusionProcess(x, f, σ, @set(x[1] >= 0 && x[2] >= 0))
CP = ControlProcess(MP, Tf, u, t, @set(u[1] >= 0 && u[1] <= 1), nothing, nothing, LagrangeMayer(lagrange, 0), 0)

bnds, V = optimal_control(CP, x0, order, trange, Mosek.Optimizer; value_function = true)


N = 100
urange = 0:0.05:1.0
function control_policy(z,s)
    one_step_MPC = MarkovBounds.extended_inf_generator(MP, V(z,s)[1], t) + lagrange
    u = urange[findmin([one_step_MPC(z..., u, s) for u in urange])[2]]
    return u
end
drift(x,p,t) = [f[1](x[1],x[2]); f[2](x[1],x[2],control_policy(x[1:2], t)); lagrange(x[1],x[2],control_policy(x[1:2], t))]
drift_uncontrolled(x,p,t) = [f[1](x[1],x[2]); f[2](x[1],x[2], u_target); 0]
diffusion(x,p,t) = [g[1](x[1]); 0; 0]
sde = DE.SDEProblem(drift, diffusion, [x0..., 0], (0.0,Tf))
sde_uncontrolled = DE.SDEProblem(drift_uncontrolled, diffusion, [x0..., 0], (0.0,Tf))
sde_ensemble = DE.EnsembleProblem(sde)
sde_ensemble_uncontrolled = DE.EnsembleProblem(sde_uncontrolled)

trajectories = DE.solve(sde_ensemble, trajectories=N)
trajectories_uncontrolled = DE.solve(sde_ensemble_uncontrolled, trajectories=N)

fig, ax = subplots(3,1)
for sol in trajectories_uncontrolled
    ax[1].plot(sol.t, [u[1] for u in sol.u], color="blue", linewidth=0.5)
    ax[2].plot(sol.t, [u[2] for u in sol.u], color="blue", linewidth=0.5)
    ax[3].plot([u[1] for u in sol.u], [u[2] for u in sol.u], color = "blue", linewidth=0.5)
end
for sol in trajectories
    ax[1].plot(sol.t, [u[1] for u in sol.u], color="red", linewidth=0.5)
    ax[2].plot(sol.t, [u[2] for u in sol.u], color="red", linewidth=0.5)
    ax[3].plot([u[1] for u in sol.u], [u[2] for u in sol.u], color = "red", linewidth=0.5)
end
ax[1].plot([0,Tf], [x_target[1], x_target[1]], color="k", linewidth=1, linestyle="dashed")
ax[2].plot([0,Tf], [x_target[2], x_target[2]], color="k", linewidth=1, linestyle="dashed")
ax[2].set_xlabel("time [s]")
ax[1].set_ylabel(L"x_1")
ax[2].set_ylabel(L"x_2")
ax[3].set_xlabel(L"x_1")
ax[3].set_ylabel(L"x_2")

empiric_control_cost = sum(sol.u[end][3] for sol in trajectories)/N
display(fig)
println("Empiric control cost = ", empiric_control_cost)
println("Lower bound on control cost = ", bnds.bounds)
println("Degree of suboptimality = ", (empiric_control_cost-bnds.bounds)/bnds.bounds)
