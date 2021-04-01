abstract type MarkovProcess end

mutable struct JumpProcess <: MarkovProcess
    x # state
    a # propensities
    h # jumps
    χ # support
end

mutable struct ReactionProcess <: MarkovProcess
    ReactionSystem::ReactionSystem
    JumpProcess::JumpProcess
    species_to_index::Dict
    species_to_state::Dict
    state_to_species::Dict
    # constructor
    function ReactionProcess(rn::ReactionSystem)
        specs = species(rn)
        S = prodstoichmat(rn) - substoichmat(rn)
        n = length(specs)
        @polyvar(x[1:n])
        spec2idx = Dict(specs[i] => i for i in 1:n)
        spec2state = Dict(specs[i] => x[i] for i in 1:n)
        state2spec = Dict(x[i] => specs[i] for i in 1:n)
        props = reformat_reactions(reactions(rn), spec2idx, x)
        jumps = reformat_jumps(S, spec2idx, x)
        support = intersect([@set(x[i] >= 0) for i in 1:n]...)
        return new(rn,JumpProcess(x,props,jumps,support),spec2idx, spec2state, state2spec)
    end
end

mutable struct DiffusionProcess <: MarkovProcess
    x # state
    f # drift
    σ # diffusion matrix
    χ # support
end

mutable struct JumpDiffusionProcess <: MarkovProcess
    JumpProcess::JumpProcess
    DiffusionProcess::DiffusionProcess
end


inf_generator(MP::JumpProcess, p::Polynomial) = sum(MP.a[i]*(subs(p, MP.x => MP.h[i]) - p) for i in 1:length(MP.a))
inf_generator(MP::ReactionProcess, p::Polynomial) = inf_generator(MP.JumpProcess,p)
inf_generator(MP::DiffusionProcess, p::Polynomial) = MP.f'*∂(p,MP.x) + 1/2*sum(∂²(p,MP.x,MP.x) .* MP.σ)
inf_generator(MP::JumpDiffusionProcess, p::Polynomial) = inf_generator(MP.JumpProcess,p) + inf_generator(MP.DiffusionProcess,p)
extended_inf_generator(MP::MarkovProcess, p::Polynomial, t::PolyVar) = ∂(p,t) + inf_generator(MP, p)

function transform_state!(P::JumpProcess, x::AbstractVector, z::AbstractVector; iv::AbstractVector = 1:length(P.x))
    Π = [P.x[i] => z[i] for i in 1:length(P.x)]
    P.a = subs.(P.a, Π...)
    P.h = [subs.(h[iv], Π...) for h in P.h]
    P.χ = intersect([@set(subs.(p, Π...) >= 0) for p in P.χ.p]...)
    P.x = x
end

function transform_state!(P::ReactionProcess, x::AbstractVector, z::AbstractVector; iv::AbstractVector = 1:length(P.JumpProcess.x))
    Π = [P.JumpProcess.x[i] => z[i] for i in 1:length(P.JumpProcess.x)]
    P.species_to_state = Dict(spec => subs(P.species_to_state[spec], Π...) for spec in keys(P.species_to_state))
    P.state_to_species = Dict(P.species_to_state[spec] => spec for spec in keys(P.species_to_state))
    transform_state!(P.JumpProcess, x, z; iv = iv)
end


function rescale_state!(P::DiffusionProcess, x0)
    transform_state!(P, P.x, P.x .* x0)
    P.f ./= x0
    P.σ ./= x0*x0'
end

function rescale_state!(P::JumpProcess, x0)
    transform_state!(P, P.x, P.x .* x0)
    for i in 1:length(P.h)
        P.h[i] ./= x0
    end
end

function rescale_state!(P::JumpDiffusionProcess, x0)
    rescale_state!(P.JumpProcess.x, x0)
    rescale_state!(P.DiffusionProcess.x, x0)
end

function rescale_state!(P::ReactionProcess, x0::AbstractVector)
    transform_state!(P, P.JumpProcess.x, P.JumpProcess.x .* x0)
    for i in 1:length(P.JumpProcess.h)
        P.JumpProcess.h[i] ./= x0
    end
end

function partition_variables(B)
    m, n = size(B) # m < n if well-posed => m
    @assert(m <= n) # implies R is m × n
    Q, R = qr(B)
    ds = []
    i, k = 1, 1
    while k <= n && i <= m
        if abs(R[i,k]) >= 1e-12
            push!(ds, k)
            i += 1
        end
        k += 1
    end
    is = setdiff(1:n, ds)
    return is, ds, Q, R
end

function project_into_subspace!(P::MarkovProcess, B, f)
    iv, dv, Q, R = partition_variables(B)
    if !isempty(dv)
        z = Vector{Polynomial}(undef, size(B,2))
        @polyvar(x[1:length(iv)])
        z[iv] .= polynomial.(x)
        z[dv] .= polynomial.(round.(R[:,dv]\(Q'*f), digits = 12) - round.((R[:,dv]\R[:,iv]), digits = 12)*x)
        transform_state!(P, x, z; iv=iv)
    end
end

function project_into_subspace!(P::ReactionProcess, x0::AbstractVector)
    B = nullspace(stoichmat(P.ReactionSystem))'
    if !isempty(B)
        f = B*x0
        project_into_subspace!(P, B, f)
    end
end
