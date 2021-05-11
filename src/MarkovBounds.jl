module MarkovBounds

using MultivariateMoments, JuMP, Reexport
@reexport using Catalyst, MomentOpt, SemialgebraicSets, DynamicPolynomials

const MM = MultivariateMoments
const MO = MomentOpt

import Base: show
import LinearAlgebra: qr, nullspace, diag, Diagonal
import Parameters: @unpack



export MarkovProcess, JumpProcess, ReactionProcess, DiffusionProcess,
       JumpDiffusionProcess, ControlProcess,
       Bounds,
       ExitProbability, TerminalSetProbability, LagrangeMayer,
       extended_inf_generator, inf_generator,
       stationary_gmp, transient_gmp,
       stationary_mean, transient_mean,
       stationary_variance, transient_variance,
       stationary_probability_mass,
       optimal_control


include("processes.jl")
include("gmp.jl")
include("utils.jl")
include("distributed.jl")
end
