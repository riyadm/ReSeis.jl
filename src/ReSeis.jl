module ReSeis

using Distributions
using DataFrames
using FFTW
using RecipesBase

import StatsBase: sample     

# patterns
include("patterns.jl")
include("patternset.jl")

# simulators
include("simulators.jl")

# solvers
include("solvers.jl")

# problems
include("problems.jl")

# seismic utilities
include("seis_utils.jl")

# plot recipes
include("plotrecipes/empirical_pattern.jl")

export 
    # pattern set
    patternset,
    sample,
    
    # empirical pattern
    reals,
    realcount,
    labels,
    stack,

    # seismic related
    Kennet,
    ricker,
    d2t    
end