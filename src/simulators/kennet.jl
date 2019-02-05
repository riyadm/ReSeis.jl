struct Kennet <: AbstractSimulator
    wavelet::Vector{Float64}
    dt::Float64
    mopt::Int
    fs::Int   
end

function Kennet(;wavelet::Vector{Float64}=[], dt::Float64=1e-3, mopt::Int=2, fs::Int=0)
    @assert !isempty(wavelet) "No wavelet provided."
    @assert dt > 0 "Sampling rate dt must be positive."
    @assert mopt >= 0 "Multiple order mopt must be non-negative."
    
    Kennet(wavelet, dt, mopt, fs)
end

