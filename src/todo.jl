# todo:
#   - variance bounds
#   - det Cov
#   - probability of semialgebraic set
#   - anything else cool?
transient_gmp(MP::ReactionProcess, order::Int, x0::AbstractVector, trange::AbstractVector, solver) = transient_gmp(MP.JumpProcess, order, x0, trange, solver)

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



function build_approximate_model!(model::GMPModel, ::AbstractPrimalMode)
    degree = MomentOpt.approximation_info(model).degree
    # initiate moments
    pvar = Dict(i => MomentOpt.moments_variable(approximation_model(model), v, degree) for (i,v) in MomentOpt.gmp_variables(model))
    # add substitutions

    # add measure condition on moments
    just = Dict{Int, Any}()
    scheme_parts = Dict{Int, Vector{MomentOpt.SchemePart}}()
    for (i, v) in MomentOpt.gmp_variables(model)
        just[i] = []
        scheme_parts[i] = MomentOpt.approximation_scheme(model, v)
        for sp in scheme_parts[i]
            cref = MomentOpt.primal_scheme_constraint(approximation_model(model), sp, pvar[i])
            push!(just[i], cref)
        end
    end
    # add constraints
    pcon = Dict{Int, Vector{ConstraintRef}}()
    for (i, con) in MomentOpt.gmp_constraints(model)
        if shape(con) isa MomentOpt.MomentConstraintShape
            cref = @constraint approximation_model(model) integrate(pvar, jump_function(con)) in moi_set(con)
            pcon[i] = [cref]
        elseif shape(con) isa MeasureConstraintShape
            # add measure constraints
            refmeas = moi_set(con)
            mons = monomials(maxdegree_basis(approx_basis(refmeas), variables(refmeas), approximation_degree(model)))
            pcon[i] = @constraint approximation_model(model) sum(c.*(MM.expectation.(mons, pvar[index(v)])) for (c, v) in jump_function(con)) .== integrate.(mons, refmeas)
        end
    end
end
