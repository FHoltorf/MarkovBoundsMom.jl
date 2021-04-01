function stationary_gmp(MP::MarkovProcess, order::Int, solver)
    test_fxns = polynomial.(monomials(MP.x, 0:order))
    gmp = GMPModel(solver)
    @variable(gmp, P∞, Meas(MP.x, support = MP.χ))
    @constraint(gmp, dynamics[ϕ in test_fxns[1:end-1]], Mom(inf_generator(MP, ϕ), P∞) == 0.0)
    @constraint(gmp, normalization, Mom(test_fxns[end], P∞) == 1.0)
    return gmp
end

stationary_gmp(MP::ReactionProcess, order::Int, solver) = stationary_gmp(MP.JumpProcess, order, solver)

function transient_gmp(MP::MarkovProcess, order::Int, x0::AbstractVector, trange::AbstractVector, solver)
    nt = length(trange)
    Δt = [trange[1], [trange[i] - trange[i-1] for i in 2:nt]...]
    @polyvar(t)
    z = [MP.x..., t]
    test_fxns = polynomial.(monomials(z, 0:order))
    gmp = GMPModel(solver)
    @variable(gmp, Pₜ[i in 1:nt], Meas(MP.x, support = MP.χ))
    @variable(gmp, ξₜ[i in 1:nt], Meas(z, support = intersect(MP.χ, @set(t >= 0), @set(Δt[i] - t >= 0)) ))
    @constraint(gmp, dynamics[ϕ in test_fxns, i in 1:nt],
                     louiville(MP,Pₜ,ξₜ,ϕ,MP.x,t,i,Δt) == (i == 1 ? split_poly(ϕ,MP.x,0)(x0...) : 0.0))
    return gmp
end

transient_gmp(MP::ReactionProcess, order::Int, x0::AbstractVector, trange::AbstractVector, solver) = transient_gmp(MP.JumpProcess, order, x0, trange, solver)

function mean(gmp, x)
    P = gmp[:P∞]
    @objective(gmp, Min, Mom(polynomial(x), P))
    optimize!(gmp)
    lb = dual_objective_value(gmp)
    @objective(gmp, Max, Mom(polynomial(x), P))
    optimize!(gmp)
    ub = dual_objective_value(gmp)
    return [lb, ub]
end

function mean(gmp, x, t)
    P = gmp[:Pₜ][t]
    @objective(gmp, Min, Mom(polynomial(x), P))
    optimize!(gmp)
    lb = dual_objective_value(gmp)
    @objective(gmp, Max, Mom(polynomial(x), P))
    optimize!(gmp)
    ub = dual_objective_value(gmp)
    return [lb, ub]
end

function prepare_gmp(rn, x0, solver, scales, auto_scaling)
    P = ReactionProcess(rn)
    project_into_subspace!(P, [x0[s] for s in species(rn)])
    if auto_scaling
        scales = stoich_bounds(rn, x0, solver)
    end
    x_scale = [scales[P.state_to_species[x]] for x in P.JumpProcess.x]
    rescale_state!(P, x_scale)
    #x_init = [x0[P.state_to_species[P.JumpProcess.x]] for x in P.JumpProcess.x] ./ x_scale
    return P
end

function mean(rn::ReactionSystem, specs::AbstractVector, x0::Dict, orders::AbstractVector, solver; scales = Dict(spec => 1.0 for spec in species(rn)), auto_scaling = true)
    P = prepare_gmp(rn, x0, solver, scales, auto_scaling)
    sol = Dict()
    for m in orders
        sol[m] = []
        for s in specs
            bnds = mean(stationary_gmp(P,m,solver), P.species_to_state[s])
            push!(sol[m], bnds)
        end
    end
    return sol
end

function stoich_bounds(rn::ReactionSystem, x0::Dict, solver)
    S = stoichmat(rn)
    B = nullspace(S)
    specs = species(rn)
    x_init = [x0[spec] for spec in specs]
    m = Model(solver)
    @variable(m, x[1:length(specs)] >= 0)
    @constraint(m, B'*x .== B'*x_init)
    scales = Dict()
    for i in 1:length(specs)
        @objective(m, Max, x[i])
        optimize!(m)
        if termination_status(m) == MOI.DUAL_INFEASIBLE
            @warn  "state space unbounded, autoscaling failed!"
            scales[specs[i]] = 1.0
        else
            scales[specs[i]] = min(1.0, objective_value(m))
        end
    end
    return scales
end

# need to do
function var(gmp, x) # need workaround
    P = gmp[:P∞]
    s = try @variable(gmp, s) catch; gmp[:s] end
    @constraint(gmp, [Mom(polynomial(x^2), P)-s   Mom(polynomial(x), P)
                      Mom(polynomial(x), P)       1.0] in PSDCone())
    @objective(gmp, Min, s)
    optimize!(gmp)
    return dual_objective_value(gmp)
end

function var(gmp, x, t) # need workaround
    P = gmp[:Pₜ][t]
    s = try @variable(gmp, s) catch; gmp[:s] end
    @constraint(gmp, [Mom(polynomial(x^2), P)-s   Mom(polynomial(x), P)
                      Mom(polynomial(x), P)       1.0] in PSDCone())
    @objective(gmp, Min, s)
    optimize!(gmp)
    return dual_objective_value(gmp)
end
