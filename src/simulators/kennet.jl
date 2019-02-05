struct Kennet <: AbstractSimulator
    wavelet::Vector{Float64}
    dt::Float64
    mopt::Int
    fs::Int
    varnames::Vector{Symbol}
end

function Kennet(;wavelet::Vector{Float64}=[], dt::Float64=1e-3, 
        mopt::Int=2, fs::Int=0, varnames::Vector{Symbol}=[:rho, :vp, :dh])
    @assert !isempty(wavelet) "No wavelet provided."
    @assert dt > 0 "Sampling rate dt must be positive."
    @assert mopt >= 0 "Multiple order mopt must be non-negative."
    
    Kennet(wavelet, dt, mopt, fs, varnames)
end

(sim::Kennet)(lyr::Dict{Symbol,Vector{Float64}}) = simkennet(lyr, wavelet=s.wavelet, dt=s.dt, mopt=s.mopt, fs=s.fs, varnames=s.varnames)
(sim::Kennet)(p::AbstractPattern) = sim(sample(p))
(sim::Kennet)(ps::Vector{<:AbstractPattern}) = sim.(ps)
(sim::Kennet)(pset::PatternSet) = sim(sample(pset))