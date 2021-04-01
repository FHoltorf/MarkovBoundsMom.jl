stoichmat(rn::ReactionSystem) = prodstoichmat(rn) - substoichmat(rn)

reformat_jumps(S::Matrix, species_to_index::Dict, x::AbstractVector) =  [x .+ S[i,:] for i in 1:size(S,1)]

split_poly(p::Polynomial, x::Vector{PolyVar}, t) = prod(x .^ p.x.Z[1][1:length(x)])*t^p.x.Z[1][end]

Base.show(io::IO, p::MarkovProcess) = println(io, "MarkovProcess")

∂(p,x) = differentiate.(p,x)
∂²(p,x,y) = differentiate(∂(p,x),y)

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
