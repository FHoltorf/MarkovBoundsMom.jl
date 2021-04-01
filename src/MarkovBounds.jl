module MarkovBounds

using Catalyst, DynamicPolynomials, MomentOpt, JuMP
import Base: show
import LinearAlgebra: qr, nullspace
import Parameters: @unpack

export MarkovProcess, JumpProcess, ReactionProcess, DiffusionProcess, JumpDiffusionProcess,
       stationary_gmp, transient_gmp, mean, var

include("processes.jl")
include("gmp.jl")
include("utils.jl")
end
