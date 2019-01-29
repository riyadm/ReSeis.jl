module ReSeis

using Distributions
using DataFrames
using GaussianMixtures
using FFTW

import StatsBase: sample
import Base: show      

include("./seis_utils.jl")
include("./patterns.jl")
include("./patternset.jl")

export 
    patternset,
    sample,
    getreals,
    realcount,
    labels,
    ricker,
    synthetic

end