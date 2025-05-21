"""
Main module for `SideKicks.jl` -- an analysis package for performing parameter inference on stellar binaries containing compact objects.

"""
module SideKicks

include("Constants.jl")
include("Orbits.jl")
include("BHFallback.jl")
include("Observations.jl")
include("Priors.jl")
include("TuringModels.jl")
include("KickMCMC.jl")
#include("CornerPlots.jl")
include("OutputStorage.jl")

end # module
