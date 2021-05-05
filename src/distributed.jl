using LightGraphs, MetaGraphs, SumOfSquares

export grid_graph, Partition

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
    n_eff = n .- inf_top .- inf_floor
    Δx = zeros(length(n))
    for k in 1:length(n)
        if n_eff[k] == 0
            lb[k] = (lb[k] + ub[k])/2
            ub[k] = lb[k]
            Δx[k] = 0
        else
            Δx[k] = (ub[k] - lb[k]) / n_eff[k]
        end
    end
    nv = prod(n)
    mg = MetaGraph(SimpleGraph(nv))
    for i in 1:nv
        idx = invert_index(i, n)
        subset = FullSpace()
        for k in 1:length(n)
            if n_eff[k] >= 0
                if inf_top[k] == 1 && idx[k] == n[k]
                    subset = intersect(subset, @set(x[k] >= ub[k]))
                elseif inf_floor[k] == 1 && idx[k] == 1
                    subset = intersect(subset, @set(x[k] <= lb[k]))
                else
                    subset = intersect(subset, @set(x[k] >= lb[k] + (idx[k]-1)*Δx[k] && x[k] <= lb[k] + idx[k]*Δx[k]))
                end
            end
        end
        set_prop!(mg, i, :cell, subset)
        for k in 1:length(n)
            if idx[k] < n[k] && n_eff[k] >= 0
                idx[k] += 1
                j = linearize_index(idx, n)
                add_edge!(mg, i, j)
                #set_prop!(mg, Edge(i,j), :overlap, get_face(subset, - x[k] + (lb[k] + (idx[k]-1)*Δx[k]), [x[k] - (lb[k] + (idx[k] - 2)*Δx[k])]))
                set_prop!(mg, Edge(i,j), :interface, (k, lb[k] + (idx[k]-1)*Δx[k]))
                idx[k] -= 1
            end
        end
    end
    return mg, x -> linearize_index(min.(floor.(Int64, x ./Δx) .+ 1, n), n)
end

#=
function get_face(subset, poly, obs_polys = [])
    return basicsemialgebraicset(algebraicset([poly]), [p for p in subset.p if !(p in obs_polys) && p != poly])
end
=#

function control_gmp(CP::ControlProcess, x0::AbstractVector, order::Int, trange::AbstractVector, p::Partition, solver)
    @assert typeof(CP.Objective) == LagrangeMayer

    MP = CP.MP
    nₜ = length(trange)
    Δt = [trange[1], [trange[i] - trange[i-1] for i in 2:nₜ]...]

    model = Model(solver)

    mons_xt = monomials([MP.x...,CP.t], 0:order)
    @variable(model, w[vertices(p.graph), k in 1:nₜ], Poly(mons_xt))
    for v in vertices(p.graph), k in 1:nₜ
        lhs = extended_inf_generator(MP, w[v,k], CP.t) + CP.Objective.l
        X = props(p.graph, v)[:cell]
        coeffs = cat(X.p, [CP.t, Δt[k]-CP.t], CP.U.p, dims = 1)
        d_lhs = maxdegree(lhs)
        mons = monomials([MP.x...,CP.u...,CP.t], 0:ceil(Int64, d_lhs/2))
        rhs = @variable(model, [1], SOSPoly(mons))[1]
        for c in coeffs
            d_rhs = ceil(Int64, (d_lhs - maxdegree(c))/2)
            mons = monomials([MP.x...,CP.u...,CP.t], 0:d_rhs)
            sos = @variable(model, [1], SOSPoly(mons))[1]
            rhs += c*sos
        end
        @constraint(model, lhs == rhs)
    end

    for v in vertices(p.graph), k in 2:nₜ
        lhs = subs(w[v,k], CP.t => 0) - subs(w[v,k-1], CP.t => Δt[k])
        d_lhs = maxdegree(lhs)
        mons = monomials(MP.x, 0:ceil(Int64, d_lhs/2))
        rhs = @variable(model, [1], SOSPoly(mons))[1]
        X = props(p.graph, v)[:cell]
        for c in X.p
            d_rhs = ceil(Int64, (d_lhs - maxdegree(c))/2)
            mons = monomials(MP.x, 0:d_rhs)
            sos = @variable(model, [1], SOSPoly(mons))[1]
            rhs += c*sos
        end
        @constraint(model, lhs == rhs)
    end

    for e in edges(p.graph), k in 1:nₜ
        v, val = props(p.graph, e)[:interface]
        @constraint(model, subs(w[e.src,k] - w[e.dst,k], MP.x[v] => val) == 0)
    end

    for v in vertices(p.graph)
        lhs = subs(w[v,nₜ], CP.t => Δt[nₜ])
        X = props(p.graph, v)[:cell]
        d_lhs = maxdegree(lhs)
        mons = monomials(MP.x, 0:ceil(Int64, d_lhs/2))
        rhs = @variable(model, [1], SOSPoly(mons))[1]
        for c in X.p
            d_rhs = ceil(Int64, (d_lhs - maxdegree(c))/2)
            mons = monomials(MP.x, 0:d_rhs)
            sos = @variable(model, [1], SOSPoly(mons))[1]
            rhs += c*sos
        end
        @constraint(model, CP.Objective.m - lhs == rhs)
    end
    @objective(model, Max, w[p.get_vertex(x0),1](x0..., 0))
    return model
end

function optimal_control(CP::ControlProcess, x0::AbstractVector, order::Int, trange::AbstractVector, p::Partition, solver; value_function = false)
    gmp = control_gmp(CP, x0, order, trange, p, solver)
    optimize!(gmp)
    lb, stat, time = extract_solution(gmp)
    if value_function
        return Bounds(:control_problem, order, lb, InfoData(stat, time)), value_function_approximation(gmp, p, CP.t, trange)
    else
        return Bounds(:control_problem, order, lb, InfoData(stat, time))
    end
end

function value_function_approximation(gmp, p::Partition, t, trange)
    V_pieces = Dict((k, i) => subs(value(gmp[:w][k,i]), t => t - (i > 1 ? t - trange[i] : 0)) for k in vertices(p.graph), i in 1:length(trange))
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
