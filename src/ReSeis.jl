module ReSeis

using Distributions
using DataFrames
using FFTW
using GaussianMixtures
using RecipesBase

import StatsBase: sample     

# patterns
include("patterns.jl")
include("patternset.jl")

# seismic utilities
include("seis_utils.jl")

# plot recipes
include("plotrecipes/patterns.jl")

export 
    patternset,
    sample,
    getreals,
    realcount,
    labels,

    # seismic related
    ricker,
    synthetic,
    sampleresponse,
    d2t
end