stoichmat(rn::ReactionSystem) = prodstoichmat(rn) - substoichmat(rn)

reformat_jumps(S::Matrix, species_to_index::Dict, x::AbstractVector) =  [x .+ S[i,:] for i in 1:size(S,1)]

split_poly(p::Polynomial, x::AbstractVector, t) = prod(x .^ p.x.Z[1][1:length(x)])*t^p.x.Z[1][end]

Base.show(io::IO, p::MarkovProcess) = println(io, "MarkovProcess")

∂(p,x) = differentiate.(p,x)
∂²(p,x,y) = differentiate(∂(p,x),y)
∂(X) = [intersect(@set(X.p[i] == 0), [@set(X.p[k] >= 0) for k in 1:length(X.p) if k != i]...) for i in 1:length(X.p)]

function reformat_reactions(rxns::Vector{Reaction}, species_to_index::Dict, x::AbstractVector)
    props = []
    for r in rxns
        @unpack rate, substrates, substoich, only_use_rate = r
        a = rate*polynomial(MonomialVector(x,0))
        if !only_use_rate
            for (s, ν) in enumerate(substoich)
                idx = species_to_index[substrates[s]]
                a *= prod(x[idx] - i for i in 0:ν-1)/factorial(ν) # consistent with rxns
            end
        end
        push!(props, a)
    end
    return props
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

extract_primal_solution(gmp::GMPModel) = objective_value(gmp), termination_status(gmp), MOI.get(gmp.approximation_model, MOI.SolveTime())
extract_solution(gmp::GMPModel) = dual_objective_value(gmp), dual_status(gmp), MOI.get(gmp.approximation_model, MOI.SolveTime())
extract_primal_solution(model::Model) = objective_value(model), termination_status(model), MOI.get(model, MOI.SolveTime())
extract_solution(model::Model) = dual_objective_value(model), dual_status(model), MOI.get(model, MOI.SolveTime())

function remap_states!(P::ReactionProcess, bnds::Bounds)
    bnds.bounds = Dict(P.state_to_species[key] => bnds.bounds[key] for key in keys(bnds.bounds))
    bnds.info.solution_time = Dict(P.state_to_species[key] => bnds.info.solution_time[key] for key in keys(bnds.info.solution_time))
    bnds.info.termination_status = Dict(P.state_to_species[key] => bnds.info.termination_status[key] for key in keys(bnds.info.termination_status))
end

function build_approximate_model(model::GMPModel, ::MO.AbstractPrimalMode)
    degree = MO.approximation_info(model).degree
    amodel = Model()
    # initiate moments
    pvar = Dict(i => MO.moments_variable(amodel, v, degree) for (i,v) in MO.gmp_variables(model))
    # add substitutions

    # add measure condition on moments
    just = Dict{Int, Any}()
    scheme_parts = Dict{Int, Vector{MO.SchemePart}}()
    for (i, v) in MO.gmp_variables(model)
        just[i] = []
        scheme_parts[i] = MO.approximation_scheme(model, v)
        for sp in scheme_parts[i]
            cref = MO.primal_scheme_constraint(amodel, sp, pvar[i])
            push!(just[i], cref)
        end
    end

    # add constraints
    pcon = Dict{Int, Vector{ConstraintRef}}()
    for (i, con) in MO.gmp_constraints(model)
        if shape(con) isa MO.MomentConstraintShape
            cref = @constraint amodel integrate(pvar, jump_function(con)) in moi_set(con)
            pcon[i] = [cref]
        elseif shape(con) isa MeasureConstraintShape
            # add measure constraints
            refmeas = moi_set(con)
            mons = monomials(maxdegree_basis(approx_basis(refmeas), variables(refmeas), approximation_degree(model)))
            pcon[i] = @constraint amodel sum(c.*(MultivariateMoments.expectation.(mons, pvar[index(v)])) for (c, v) in jump_function(con)) .== integrate.(mons, refmeas)
        end
    end
    return amodel, pvar
end

function get_moment(p::Polynomial, meas, model)
    expr = AffExpr(0.0)
    for i in 1:length(p.a)
        index = findfirst(isequal(p.x[i]), meas.x)
        add_to_expression!(expr, p.a[i], VariableRef(model, MOI.VariableIndex(index)))
    end
    return expr
end

function UpperDiag(vars, n)
    M = AffExpr.(zeros(n,n))
    k = 1
    for i in 1:n, j in i:n
        M[i,j] += vars[k]
        k += 1
    end
    return M
end
