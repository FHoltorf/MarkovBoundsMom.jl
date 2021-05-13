using LightGraphs, MetaGraphs, SumOfSquares

export grid_graph, Partition, props

struct Partition
    graph
    get_vertex
end

linearize_index(idx,rs) = idx[1] + sum((idx[i] - 1) * prod(rs[1:i-1]) for i in 2:length(idx))

function invert_index(idx, rs)
    inv_idx = similar(rs)
    for i in length(rs):-1:2
        fac = prod(rs[1:i-1])
        n = div(idx - 1, fac)
        inv_idx[i] = n + 1
        idx -= n*fac
    end
    inv_idx[1] = idx
    return inv_idx
end

function grid_graph(x, lb, ub, n; inf_top = zeros(Int64, length(ub)), inf_floor = zeros(Int64, length(lb)))
    degen_dim = []
    nondegen_dim = []
    for k in 1:length(n)
        if n[k] == 1
            push!(degen_dim, k)
        else
            push!(nondegen_dim, k)
        end
    end
    n_eff = n .- inf_top .- inf_floor
    Δx = zeros(length(n))
    for k in nondegen_dim
        if n_eff[k] == 0
            lb[k] = (lb[k] + ub[k])/2
            ub[k] = lb[k]
            Δx[k] = 1.0
        else
            Δx[k] = (ub[k] - lb[k]) / n_eff[k]
        end
    end
    nv = prod(n)
    mg = MetaGraph(SimpleGraph(nv))
    for i in 1:nv
        idx = invert_index(i, n)
        subset = FullSpace()
        for k in nondegen_dim
            if inf_top[k] == 1 && idx[k] == n[k]
                subset = intersect(subset, @set(x[k] >= ub[k]))
            elseif inf_floor[k] == 1 && idx[k] == 1
                subset = intersect(subset, @set(x[k] <= lb[k]))
            else
                subset = intersect(subset, @set(x[k] >= lb[k] + (idx[k]-inf_floor[k]-1)*Δx[k] && x[k] <= lb[k] + (idx[k]-inf_floor[k])*Δx[k]))
            end
        end
        for k in degen_dim
            if inf_top[k] == 0
                subset = intersect(subset, @set(x[k] <= ub[k]))
            end
            if inf_floor[k] == 0
                subset = intersect(subset, @set(x[k] >= lb[k]))
            end
        end
        set_prop!(mg, i, :cell, subset)
        for k in nondegen_dim
            if idx[k] < n[k]
                idx[k] += 1
                j = linearize_index(idx, n)
                add_edge!(mg, i, j)
                set_prop!(mg, Edge(i,j), :interface, (k, lb[k] + (idx[k]-1-inf_floor[k])*Δx[k]))
                idx[k] -= 1
            end
        end
    end
    get_vertex = function (x)
                    idx = ones(Int64,length(n))
                    for k in nondegen_dim
                        if inf_floor[k] == 1
                            idx[k] = x[k] <= lb[k] ? 1 : min(ceil(Int64, (x[k] - lb[k])/Δx[k]) + 1, n[k])
                        else
                            idx[k] = min(floor(Int64, (x[k] - lb[k])/Δx[k]) + 1, n[k])
                        end
                    end
                    return linearize_index(idx, n)
                 end
    return mg, get_vertex
end

function grid_graph(x, x_ranges)
    n = length.(x_ranges) .- 1
    nv = prod(n)

    mg = MetaGraph(SimpleGraph(nv))
    for i in 1:nv
        idx = invert_index(i, n)
        subset = FullSpace()
        for k in 1:length(idx)
            if !isinf(x_ranges[k][idx[k]])
                subset = intersect(subset, @set(x[k] >= x_ranges[k][idx[k]]))
            end
            if !isinf(x_ranges[k][idx[k]+1])
                subset = intersect(subset, @set(x[k] <= x_ranges[k][idx[k]+1]))
            end
        end
        set_prop!(mg, i, :cell, subset)
        for k in 1:length(idx)
            if idx[k] < n[k]
                idx[k] += 1
                j = linearize_index(idx, n)
                add_edge!(mg, i, j)
                set_prop!(mg, Edge(i,j), :interface, (k, x_ranges[k][idx[k]]))
                idx[k] -= 1
            end
        end
    end
    get_vertex = function (x)
                    idx = zeros(Int64, length(n))
                    for k in 1:length(n)
                        j = findfirst(m -> m >= x[k], x_ranges[k]) - 1
                        if isnothing(j)
                            @error "state outside domain"
                        end
                        idx[k] = max(j, 1)
                    end
                    return linearize_index(idx, n)
                 end
    return mg, get_vertex
end

function control_gmp(CP::ControlProcess, x0::AbstractVector, order::Int, trange::AbstractVector, p::Partition, solver)
    @assert typeof(CP.Objective) == LagrangeMayer

    MP = CP.MP
    nₜ = length(trange)
    Δt = [trange[1], [trange[i] - trange[i-1] for i in 2:nₜ]...]

    model = SOSModel(solver)

    @variable(model, w[vertices(p.graph), k in 1:nₜ], Poly(monomials([MP.x..., CP.t], 0:order)))

    for v in vertices(p.graph), k in 1:nₜ
        X = intersect(props(p.graph, v)[:cell], @set(CP.t >= 0 && 1-CP.t >= 0), CP.U)
        @constraint(model, extended_inf_generator(MP, w[v,k], CP.t; scale = Δt[k]) + Δt[k]*CP.Objective.l >= 0, domain = X)
    end

    for v in vertices(p.graph), k in 2:nₜ
        @constraint(model, subs(w[v,k], CP.t => 0) - subs(w[v,k-1], CP.t => 1) >= 0, domain=props(p.graph, v)[:cell])
    end

    for e in edges(p.graph), k in 1:nₜ
        i, val = props(p.graph, e)[:interface]
        if all(subs(MP.σ, MP.x[i] => val) .== 0)
            X = FullSpace()
            for p in props(p.graph, e.dst)[:cell].p
                if !(variables(p) == [MP.x[i]])
                    p = subs(p, MP.x[i] => val)
                    X = intersect(X, @set(p >= 0))
                end
            end
            X = intersect(X, CP.U, @set(CP.t >= 0 && 1-CP.t >= 0))
            @constraint(model, subs((w[e.dst,k] - w[e.src,k])*MP.f[i], MP.x[i] => val) >= 0, domain=X)
        else
            @constraint(model, subs(w[e.src,k] - w[e.dst,k], MP.x[i] => val) == 0)
        end
    end
    for v in vertices(p.graph)
        @constraint(model, CP.Objective.m - subs(w[v,nₜ], CP.t => 1) >= 0, domain = props(p.graph, v)[:cell])
    end
    @objective(model, Max, w[p.get_vertex(x0),1](x0..., 0))
    return model
end

function optimal_control(CP::ControlProcess, x0::AbstractVector, order::Int, trange::AbstractVector, p::Partition, solver; value_function = false)
    gmp = control_gmp(CP, x0, order, trange, p, solver)
    optimize!(gmp)
    lb, stat, time = extract_primal_solution(gmp)
    if value_function
        return Bounds(:control_problem, order, lb, InfoData(stat, time)), value_function_approximation(gmp, p, CP.t, trange)
    else
        return Bounds(:control_problem, order, lb, InfoData(stat, time))
    end
end

function value_function_approximation(gmp, p::Partition, t, trange)
    Δt = [trange[1], [trange[i] - trange[i-1] for i in 2:length(trange)]...]
    V_pieces = Dict((k, i) => subs(value(gmp[:w][k,i]), t => (t - (i > 1 ? trange[i] : 0))/Δt[i]) for k in vertices(p.graph), i in 1:length(trange))
    V = function (x,t)
            if t > trange[end] || t < 0
                @warn string("t = ", t, " outside the domain of V")
                i = t < 0 ? 0 : length(trange) - 1
            else
                i = findlast(ti -> ti <= t, trange[1:end-1])
            end
            k = p.get_vertex(x)
            return isnothing(i) ? (V_pieces[k,1], V_pieces[k,1](x...,t)) : (V_pieces[k,i+1], V_pieces[k,i+1](x...,t))
        end
    return V
end

function inf_horizon_control_gmp(CP::ControlProcess, x0::AbstractVector, order::Int, trange::AbstractVector, p::Partition, solver; value_function = false)
    @assert typeof(CP.Objective) == LagrangeMayer
    @assert CP.Objective.m == 0

    MP = CP.MP
    nₜ = length(trange)
    Δt = [trange[1], [trange[i] - trange[i-1] for i in 2:nₜ]...]
    ρ = CP.discount_factor
    model = SOSModel(solver)

    @variable(model, w[vertices(p.graph), k in 1:nₜ], Poly(monomials([MP.x..., CP.t], 0:order)))
    @variable(model, w∞[vertices(p.graph)], Poly(monomials(MP.x, 0:order)))

    for v in vertices(p.graph), k in 1:nₜ
        X = intersect(props(p.graph, v)[:cell], @set(CP.t >= 0 && 1-CP.t >= 0), CP.U)
        @constraint(model, extended_inf_generator(MP, w[v,k], CP.t; scale = Δt[k]) - 2*ρ*w[v,k] + Δt[k]*CP.Objective.l >= 0, domain = X)
    end

    for v in vertices(p.graph), k in 2:nₜ
        @constraint(model, subs(w[v,k], CP.t => 0) - subs(w[v,k-1], CP.t => 1) >= 0, domain=props(p.graph, v)[:cell])
    end

    for e in edges(p.graph), k in 1:nₜ
        i, val = props(p.graph, e)[:interface]
        if all(subs(MP.σ, MP.x[i] => val) .== 0)
            X = FullSpace()
            for p in props(p.graph, e.dst)[:cell].p
                if !(variables(p) == [MP.x[i]])
                    p = subs(p, MP.x[i] => val)
                    X = intersect(X, @set(p >= 0))
                end
            end
            X = intersect(X, CP.U, @set(CP.t >= 0 && 1-CP.t >= 0))
            @constraint(model, subs((w[e.dst,k]-w[e.src,k])*MP.f[i], MP.x[i] => val) >= 0, domain=X)
        else
            @constraint(model, subs(w[e.src,k]-w[e.dst,k], MP.x[i] => val) == 0)
        end
    end

    for v in vertices(p.graph)
        X = intersect(props(p.graph, v)[:cell], @set(CP.t >= 0), CP.U)
        @constraint(model, extended_inf_generator(MP, w∞[v], CP.t) - 2*ρ*w∞[v] + CP.Objective.l >= 0, domain = X)
        @constraint(model, w∞[v] - subs(w[v,nₜ], CP.t => 1) >= 0, domain=props(p.graph, v)[:cell])
    end

    for e in edges(p.graph), k in 1:nₜ
        i, val = props(p.graph, e)[:interface]
        if all(subs(MP.σ, MP.x[i] => val) .== 0)
            X = FullSpace()
            for p in props(p.graph, e.dst)[:cell].p
                if !(variables(p) == [MP.x[i]])
                    p = subs(p, MP.x[i] => val)
                    X = intersect(X, @set(p >= 0))
                end
            end
            X = intersect(X, CP.U, @set(CP.t >= 0))
            @constraint(model, subs((w∞[e.dst]-w∞[e.src])*MP.f[i], MP.x[i] => val) >= 0, domain=X)
        else
            @constraint(model, subs(w∞[e.src]-w∞[e.dst], MP.x[i] => val) == 0)
        end
    end
    @objective(model, Max, w[p.get_vertex(x0),1](x0..., 0))
    return model
end

function inf_horizon_control(CP::ControlProcess, x0::AbstractVector, order::Int, trange::AbstractVector, p::Partition, solver; value_function = false)
    gmp = inf_horizon_control_gmp(CP, x0, order, trange, p, solver)
    optimize!(gmp)
    lb, stat, time = extract_primal_solution(gmp)
    if value_function
        return Bounds(:control_problem, order, lb, InfoData(stat, time)), value_function_approximation_inf_horizon(gmp, p, CP.t, trange)
    else
        return Bounds(:control_problem, order, lb, InfoData(stat, time))
    end
end

function value_function_approximation_inf_horizon(gmp, p::Partition, t, trange)
    nt = length(trange)
    Δt = [trange[1], [trange[i] - trange[i-1] for i in 2:nt]...]
    V_pieces = Dict((k, i) => subs(value(gmp[:w][k,i]), t => (t - (i > 1 ? trange[i] : 0))/Δt[i]) for k in vertices(p.graph), i in 1:length(trange))
    merge!(V_pieces, Dict((k, nt+1) => subs(value(gmp[:w∞][k]), t => t - trange[end])+0.0*t for k in vertices(p.graph)))
    V = function (x,t)
            k = p.get_vertex(x)
            if t < 0
                @warn string("t = ", t, " outside the domain of V")
                i = t < 0 ? 0 : length(trange) - 1
            else
                i = findlast(ti -> ti <= t, trange)
            end
            return isnothing(i) ? (V_pieces[k,1], V_pieces[k,1](x...,t)) : (V_pieces[k,i+1], V_pieces[k,i+1](x...,t))
        end
    return V
end
