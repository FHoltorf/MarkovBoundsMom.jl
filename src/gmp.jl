## Stationary problems
function stationary_gmp(MP::MarkovProcess, order::Int, solver)
    test_fxns = polynomial.(monomials(MP.x, 1:order))
    gmp = GMPModel(solver)
    set_approximation_mode(gmp, PRIMAL_RELAXATION_MODE())
    @variable(gmp, P∞, Meas(MP.x, support = MP.X))
    @constraint(gmp, dynamics[ϕ in test_fxns], Mom(inf_generator(MP, ϕ), P∞) == 0.0)
    @constraint(gmp, normalization, Mom(1, P∞) == 1.0)
    return gmp
end

stationary_gmp(MP::ReactionProcess, order::Int, solver) = stationary_gmp(MP.JumpProcess, order, solver)

function stationary_moment(gmp, x)
    P = gmp[:P∞]
    @objective(gmp, Min, Mom(polynomial(x), P))
    optimize!(gmp)
    lb, stat_lb, time_lb = extract_solution(gmp)
    @objective(gmp, Max, Mom(polynomial(x), P))
    optimize!(gmp)
    ub, stat_ub, time_ub = extract_solution(gmp)
    return [lb, ub], [stat_lb, stat_ub], [time_lb, time_ub]
end

function prepare_reaction_system(rn::ReactionSystem, x0::Dict, solver, scales::Dict, auto_scaling)
    P = ReactionProcess(rn)
    project_into_subspace!(P, [x0[s] for s in species(rn)])
    if auto_scaling
        scales = stoich_bounds(rn, x0, solver)
    end
    x_scale = [scales[P.state_to_species[x]] for x in P.JumpProcess.x]
    x0 = [x0[P.state_to_species[x]]/scales[P.state_to_species[x]] for x in P.JumpProcess.x]
    rescale_state!(P, x_scale)
    return P, x0
end

function stationary_mean(MP::MarkovProcess, states::AbstractVector, order::Int, solver; G = [], f = [])
    if !isempty(G)
        project_into_subspace!(MP, G, f)
    end
    gmp = stationary_gmp(MP, order, solver)
    res = Dict(s => stationary_moment(gmp, s) for s in states)
    sol = Bounds((:mean, Inf), order, Dict(s => res[s][1] for s in states), InfoData(Dict(s => res[s][2] for s in states), Dict(s => res[s][3] for s in states)))
    return sol
end

function stationary_mean(rn::ReactionSystem, specs::AbstractVector, x0::Dict, order::Int, solver;
                         scales = Dict(spec => 1.0 for spec in species(rn)), auto_scaling = false)
    P, x0 = prepare_reaction_system(rn, x0, solver, scales, auto_scaling)
    sol = stationary_mean(P.JumpProcess, [P.species_to_state[s] for s in specs], order, solver)
    remap_states!(P, sol)
    return sol
end

function transient_gmp(MP::MarkovProcess, x0::AbstractVector, order::Int, trange::AbstractVector, solver)
    nt = length(trange)
    Δt = [trange[1], [trange[i] - trange[i-1] for i in 2:nt]...]
    @polyvar(t)
    z = [MP.x..., t]
    test_fxns = polynomial.(monomials(z, 0:order))
    gmp = GMPModel(solver)
    set_approximation_mode(gmp, PRIMAL_RELAXATION_MODE())
    @variable(gmp, Pₜ[i in 1:nt], Meas(MP.x, support = MP.X))
    @variable(gmp, ξₜ[i in 1:nt], Meas(z, support = intersect(MP.X, @set(t >= 0), @set(Δt[i] - t >= 0))))
    @constraint(gmp, dynamics[ϕ in test_fxns, i in 1:nt],
                     louiville(MP,Pₜ[i],ξₜ[i],ϕ,t,Δt[i]) == (i == 1 ? split_poly(ϕ,MP.x,0)(x0...) : Mom(split_poly(ϕ, MP.x, 0), Pₜ[i-1])))
    return gmp
end

transient_gmp(MP::ReactionProcess, x0::AbstractVector, order::Int, trange::AbstractVector, solver) = transient_gmp(MP.JumpProcess, x0, order, trange, solver)

function transient_moment(gmp, x, idx)
    P = gmp[:Pₜ][idx]
    @objective(gmp, Min, Mom(polynomial(x), P))
    optimize!(gmp)
    lb, stat_lb, time_lb = extract_solution(gmp)
    @objective(gmp, Max, Mom(polynomial(x), P))
    optimize!(gmp)
    ub, stat_ub, time_ub = extract_solution(gmp)
    return [lb, ub], [stat_lb, stat_ub], [time_lb, time_ub]
end

function transient_mean(MP::MarkovProcess, states::AbstractVector, x0::AbstractVector, order::Int, trange::AbstractVector, solver; G = [], f = [], t = 0)
    if !isempty(G)
        project_into_subspace!(MP, G, f)
    end
    if t != 0
        @assert t in trange "t not an element of trange"
        idx = findfirst(isequal(t), trange)
    else
        idx = length(trange)
    end
    gmp = transient_gmp(MP, x0, order, trange, solver)
    res = Dict(s => transient_moment(gmp, s, idx) for s in states)
    sol = Bounds((:mean, trange[idx]), order, Dict(s => res[s][1] for s in states), InfoData(Dict(s => res[s][2] for s in states), Dict(s => res[s][3] for s in states)))
    return sol
end

function transient_mean(rn::ReactionSystem, specs::AbstractVector, x0::Dict, order::Int, trange::AbstractVector, solver;
                        scales = Dict(spec => 1.0 for spec in species(rn)), auto_scaling = false, t=0)
    P, x0 = prepare_reaction_system(rn, x0, solver, scales, auto_scaling)
    sol = transient_mean(P.JumpProcess, [P.species_to_state[s] for s in specs], x0, order, trange, solver, t=t)
    remap_states!(P, sol)
    return sol
end

# variance
function stationary_variance(MP::MarkovProcess, states::AbstractVector, order::Int, solver; G = [], f = [])
    if !isempty(G)
        project_into_subspace!(MP, G, f)
    end
    gmp = stationary_gmp(MP, order, solver)
    model, measures = build_approximate_model(gmp, approximation_mode(gmp))
    set_optimizer(model, solver)
    P = measures[1]
    res = Dict()
    for s in states
        slack = @variable(model)
        @constraint(model, [MM.expectation(s^2,P)-slack MM.expectation(s,P);
                            MM.expectation(s,P) 1] in PSDCone())
        @objective(model, Max, slack)
        optimize!(model)
        res[s] = extract_solution(model)
    end
    sol = Bounds((:variance, Inf), order, Dict(s => res[s][1] for s in states), InfoData(Dict(s => res[s][2] for s in states), Dict(s => res[s][3] for s in states)))
    return sol
end

function stationary_variance(rn::ReactionSystem, specs::AbstractVector, x0::Dict, order::Int, solver;
                             scales = Dict(spec => 1.0 for spec in species(rn)), auto_scaling = false)
    P, x0 = prepare_reaction_system(rn, x0, solver, scales, auto_scaling)
    sol = stationary_variance(P.JumpProcess, [P.species_to_state[s] for s in specs], order, solver)
    remap_states!(P, sol)
    return sol
end

function transient_variance(MP::MarkovProcess, states::AbstractVector, x0::AbstractVector, order::Int, trange::AbstractVector, solver; G = [], f = [], t=0)
    if !isempty(G)
        project_into_subspace!(MP, G, f)
    end
    if t != 0
        @assert t in trange "t not an element of trange"
        idx = findfirst(isequal(t), trange)
    else
        idx = length(trange)
    end
    gmp = transient_gmp(MP, x0, order, trange, solver)
    model, measures = build_approximate_model(gmp, approximation_mode(gmp))
    set_optimizer(model, solver)
    P = measures[idx]
    res = Dict()
    for s in states
        slack = @variable(model)
        @constraint(model, [MM.expectation(s^2,P)-slack MM.expectation(s,P);
                            MM.expectation(s,P) 1] in PSDCone())
        @objective(model, Max, slack)
        optimize!(model)
        res[s] = extract_solution(model)
    end
    sol = Bounds((:variance, trange[idx]), order, Dict(s => res[s][1] for s in states), InfoData(Dict(s => res[s][2] for s in states), Dict(s => res[s][3] for s in states)))
    return sol
end

function transient_variance(rn::ReactionSystem, specs::AbstractVector, x0::Dict, order::Int, trange::AbstractVector, solver;
                            scales = Dict(spec => 1.0 for spec in species(rn)), auto_scaling = false, t=0)
    P, x0 = prepare_reaction_system(rn, x0, solver, scales, auto_scaling)
    sol = transient_variance(P.JumpProcess, [P.species_to_state[s] for s in specs], x0, order, trange, solver, t=t)
    remap_states!(P, sol)
    return sol
end

# Probability mass
function stationary_probability_mass(event, P::MarkovProcess, order::Int, solver)
    test_fxns = polynomial.(monomials(MP.x, 1:order))
    gmp = GMPModel(solver)
    set_approximation_mode(gmp, PRIMAL_RELAXATION_MODE())
    @variable(gmp, P∞, Meas(P.x, support = P.X))
    @variable(gmp, Q∞, Meas(P.x, support = event))
    @constraint(gmp, dynamics[ϕ in test_fxns], Mom(inf_generator(MP, ϕ), P∞ + Q∞) == 0.0)
    @constraint(gmp, normalization, Mom(1, P∞ + Q∞) == 1.0)
    @objective(gmp, Max, Mom(1, Q∞))
    optimize!(m)
    ub, stat_ub, time_ub = extract_solution(gmp)
    @objective(gmp, Min, Mom(1,ν∞))
    optimize!(m)
    lb, stat_lb, time_lb = extract_solution(gmp)
    return Bounds(:probability_mass, order, Dict(event => [lb,ub]), InfoData([stat_lb,stat_ub], [time_lb,time_ub]))
end

# volume of stationary ellipsoid
function stationary_confidence_ellipsoid(MP::MarkovProcess, states, order::Int, solver)
    # vol of confidence ellipsoide ~ det(𝔼[xxᵀ])
    gmp = stationary_gmp(MP, order, solver)
    m, measures = build_approximate_model(gmp, approximation_mode(gmp))
    set_optimizer(m, solver)
    P = measures[1]
    n = length(states)
    @variable(m, z[1:div(n*(n+1),2)])
    @variable(m, G[1:n,1:n])
    Z = UpperDiag(z, n)
    U = MM.expectation.(states*states', P)
    v = MM.expectation.(states, P)
    vᵀ = transpose(v)
    d = Diagonal(diag(Z))
    Zᵀ = transpose(Z)
    @constraint(m, Schur_Cov, [U .- G v;
                               vᵀ 1] in PSDCone())
    @constraint(m, Schur_logdet, [G  Z;
                                  Zᵀ d] in PSDCone())
    @variable(m, t[1:n])
    @constraint(m, logsum[i in 1:n], [t[i], 1, Z[i,i]] in MOI.ExponentialCone())
    @objective(m, Max, sum(t[i] for i in 1:n))
    optimize!(m)
    ub, stat, time = extract_solution(m)
    return Bounds(:volume_confidence_ellipsoid, order, Dict(states => exp(ub)), InfoData(stat, time))
end

function stationary_confidence_ellipsoid(rn::ReactionSystem, x0::Dict, species::AbstractVector, order::Int, solver;
                                         scales = Dict(spec => 1.0 for spec in species(rn)), auto_scaling = false)
    P, x0 = prepare_reaction_system(rn, x0, solver, scales, auto_scaling)
    states = [P.species_to_state[s] for s in species]
    bnds = stationary_confidence_ellipsoid(P.JumpProcess, states, order, solver)
    return Bounds(bnds.type, order, Dict(species => bnds.bounds[states]), bnds.info)
end

## Control Processes

function add_PathChanceConstraint!(gmp, cc, test_fxns, Δt, CP)
    MP = CP.MarkovProcess
    nt = length(Δt)
    ∂X = ∂(cc.X)
    Qₜ = @variable(gmp, [i in 1:nt], Meas(MP.x, support = cc.X)) # residence measure
    Rₜ = @variable(gmp, [k in 1:length(∂X), i in 1:nt],
                        Meas([MP.x..., CP.u..., CP.t], support = intersect(∂X[k], @set(CP.t >= 0), @set(CP.t <= Δt[i])))) # exit measure
    Sₜ = @variable(gmp, [i in 1:nt], Meas([MP.x..., CP.u..., CP.t], support = intersect(CP.ChanceConstraint.X, CP.U, @set(CP.t >= 0), @set(Δt[i] - CP.t >= 0))))
    @constraint(gmp, [ϕ in test_fxns, i in 1:nt],
                     louiville(MP,Qₜ[i],Sₜ[i],ϕ,CP.t,Δt[i]) + sum(Mom(ϕ, Rₜ[k,i]) for k in length(∂X)) ==
                     (i == 1 ? split_poly(ϕ,MP.x,0)(x0...) : Mom(split_poly(ϕ, MP.x, 0), Qₜ[i-1])))
    @constraint(gmp, [1], Mom(1, Qₜ[nt]) >= 1 - cc.α)
    Tₜ = @variable(gmp, [i in 1:nt], Meas(MP.x, support = MP.X)) # slack measure
    @constraint(gmp, [i in 1:nt], Tₜ[i] + Qₜ[i] == gmp[:Pₜ][i])
end

function add_TerminalChanceConstraint!(gmp, cc, nt, CP)
    MP = CP.MarkovProcess
    Q = @variable(gmp, [nt], Meas(MP.x, support = cc.X))
    S = @variable(gmp, [nt], Meas(MP.x, support = MP.X))
    @constraint(gmp, [nt], Q[nt] + S[nt] == gmp[:Pₜ][nt])
    @constraint(gmp, [nt], Mom(1, Q[nt]) >= 1-cc.α)
end

function control_gmp(CP::ControlProcess, x0::AbstractVector, order::Int, trange::AbstractVector, solver)
    MP = CP.MP
    nt = length(trange)
    Δt = [trange[1], [trange[i] - trange[i-1] for i in 2:nt]...]
    gmp = GMPModel(solver)
    set_approximation_mode(gmp, PRIMAL_RELAXATION_MODE())
    test_fxns = polynomial.(monomials([MP.x..., CP.t], 0:order))
    if typeof(CP.Objective) == LagrangeMayer ||  typeof(CP.Objective) == TerminalSetProbability
        @variable(gmp, Pₜ[i in 1:nt], Meas(MP.x, support = MP.X))
        @variable(gmp, ξₜ[i in 1:nt], Meas([MP.x..., CP.u..., CP.t], support = intersect(MP.X, CP.U, @set(CP.t >= 0), @set(Δt[i] - CP.t >= 0))))
        @constraint(gmp, dynamics[ϕ in test_fxns, i in 1:nt],
                         louiville(MP,Pₜ[i],ξₜ[i],ϕ,CP.t,Δt[i]) == (i == 1 ? split_poly(ϕ,MP.x,0)(x0...) : Mom(split_poly(ϕ, MP.x, 0), Pₜ[i-1])))

        if !isnothing(CP.PathChanceConstraints)
            for cc in CP.PathChanceConstraints
                add_PathChanceConstraint(gmp, cc, test_fxns, Δt, CP)
            end
        end

        if !isnothing(CP.TerminalChanceConstraints)
            for cc in CP.TerminalChanceConstraints
                add_TerminalChanceConstraint!(gmp, cc, nt, CP)
            end
        end

        if typeof(CP.Objective) == LagrangeMayer
            @objective(gmp, Min, sum(Mom(CP.Objective.l, gmp[:ξₜ][i]) for i in 1:nt) + Mom(CP.Objective.m, gmp[:Pₜ][end]))
        else
            Q = @variable(gmp, [nt], Meas(CP.Objective.l, support = CP.Objective.X))
            R = @variable(gmp, [nt], Meas(CP.MP, support = MP.X))
            @constraint(gmp, R[nt] + Q[nt] == Pₜ[nt])
            @objective(gmp, Max, Mom(1, Q[nt]))
        end
    elseif typeof(CP.Objective) == ExitProbability
        ∂X = ∂(CP.Objective.X)
        Qₜ = @variable(gmp, [i in 1:nt], Meas(MP.x, support = CP.Objective.X))
        Rₜ = @variable(gmp, [k in 1:length(∂X), i in 1:nt],
                            Meas([MP.x..., CP.t], support = intersect(∂X[k], @set(CP.t >= 0), @set(CP.t <= Δt[i]))))
        Sₜ = @variable(gmp, [i in 1:nt], Meas([MP.x..., CP.u..., CP.t], support = intersect(CP.Objective.X, CP.U, @set(CP.t >= 0), @set(CP.t <= Δt[i]))))
        @constraint(gmp, dynamics[ϕ in test_fxns, i in 1:nt],
                         louiville(MP,Qₜ[i],Sₜ[i],ϕ,CP.t,Δt[i]) + sum(Mom(ϕ, Rₜ[k,i]) for k in length(∂X)) ==
                         (i == 1 ? split_poly(ϕ,MP.x,0)(x0...) : Mom(split_poly(ϕ, MP.x, 0), Qₜ[i-1])))
        @objective(gmp, Max, Mom(1, Qₜ[nt]))
    end
    return gmp
end

function value_function_approximation(gmp, t, trange)
    V_pieces = Dict(i => polynomial(0.0*t) for i in 1:length(trange))
    for idx in keys(gmp[:dynamics])
        i, ϕ = idx[2], idx[1]
        t0 = (i > 1 ? trange[i-1] : 0)
        V_pieces[i] += dual(gmp[:dynamics][idx])[1] * subs(ϕ, t => (t-t0))
    end
    V = function (x,t)
            if t > trange[end] || t < 0
                @warn string("t = ", t, " outside the domain of V")
                i = t < 0 ? 0 : length(trange) - 1
            else
                i = findlast(ti -> ti <= t, trange[1:end-1])
            end
            return isnothing(i) ? (V_pieces[1], V_pieces[1](x...,t)) : (V_pieces[i+1], V_pieces[i+1](x...,t))
        end
    return V
end

function optimal_control(CP::ControlProcess, x0::AbstractVector, order::Int, trange::AbstractVector, solver; value_function = false)
    gmp = control_gmp(CP, x0, order, trange, solver)
    optimize!(gmp)
    lb, stat, time = extract_solution(gmp)
    if value_function
        return Bounds(:control_problem, order, lb, InfoData(stat, time)), value_function_approximation(gmp, CP.t, trange)
    else
        return Bounds(:control_problem, order, lb, InfoData(stat, time))
    end
end
