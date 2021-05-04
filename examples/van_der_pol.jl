using DynamicPolynomials, MosekTools, MomentOpt, DifferentialEquations, PyPlot, LaTeXStrings, MarkovBounds
DE = DifferentialEquations
@polyvar(x[1:2])
@polyvar(u[1:1])
@polyvar(t)

μ = 1.0
f = [1/μ*x[2]; μ*(x[2] - 1/3*x[2]^3 - x[1]) + x[1]*u[1]]
g = [0; 0.7*x[1]]
σ = g*g'

Tf = 5.0
x0 = [0.5; 0.5]
U = @set(u[1] <= 1 && u[1] >= -1)

X = @set(x[1]^2 + x[2]^2 <= 1)

## data
order = 6
nt = 10
trange = range(Tf/nt, stop = Tf, length=nt)

## solve
MP = DiffusionProcess(x, f, σ, FullSpace())
CP = ControlProcess(MP, Tf, u, t, U, nothing, nothing, ExitProbability(X), 0)

bnds, V = optimal_control(CP, x0, order, trange, Mosek.Optimizer; value_function = true)


N = 500
urange = -1.0:0.05:1.0
function control_policy(z,s)
    one_step_MPC = extended_inf_generator(MP, V(z,s)[1], t)
    u = urange[findmin([one_step_MPC(z..., u, s) for u in urange])[2]]
    return u
end
drift(x,p,t) = [f[1](x[2]); f[2](x[1],x[2],control_policy(x[1:2], t))]
drift_uncontrolled(x,p,t) = [f[1](x[2]); f[2](x[1],x[2],0)]
diffusion(x,p,t) = [0; g[2](x[1])]
sde = DE.SDEProblem(drift, diffusion, x0, (0.0,Tf))
sde_uncontrolled = DE.SDEProblem(drift_uncontrolled, diffusion, x0, (0.0,Tf))
sde_ensemble = DE.EnsembleProblem(sde)
sde_ensemble_uncontrolled = DE.EnsembleProblem(sde_uncontrolled)

trajectories = DE.solve(sde_ensemble, trajectories=N, EM(), dt = 0.01)
trajectories_uncontrolled = DE.solve(sde_ensemble_uncontrolled, trajectories=N, EM(), dt = 0.01)

P_uncontrolled = 0
fig, ax = subplots()
for sol in trajectories_uncontrolled
    max_idx = findfirst(m -> m[1]^2 + m[2]^2 > 1, sol.u)
    if max_idx == nothing
        global P_uncontrolled += 1
    end
    max_idx = (max_idx == nothing ? length(sol.u) : max_idx)
    ax.plot([sol.u[i][1] for i in 1:max_idx], [sol.u[i][2] for i in 1:max_idx], color="blue", linewidth=0.5)
end

P_controlled = 0
for sol in trajectories
    max_idx = findfirst(m -> m[1]^2 + m[2]^2 > 1, sol.u)
    if max_idx == nothing
        global P_controlled += 1
    end
    max_idx = (max_idx == nothing ? length(sol.u) : max_idx)
    ax.plot([sol.u[i][1] for i in 1:max_idx], [sol.u[i][2] for i in 1:max_idx], color="red", linewidth=0.5)
end
ax.plot(cos.(-π:0.01:π), sin.(-π:0.01:π), color = "black", linewidth= 2)
ax.scatter([x0[1]], [x0[2]], marker="o", color="black", s=5)
ax.set_xlabel(L"x_1")
ax.set_ylabel(L"x_2")
display(fig)

println("controlled: P(exit) = ", P_controlled/N)
println("uncontrolled: P(exit) = ", P_uncontrolled/N)
