function ricker(f₀::Float64, Δt::Float64)
    @assert f₀ > 0 "Peak frequency f₀ must be positive!"
    @assert Δt > 0 "Time interval Δt must be positive!"
    
    # signal length
    length = 2.0 / f₀ / Δt
    
    # number of points in each direction
    n = ceil(Int, length / 2)
    
    # time axis
    t = Δt * (-n:n)
    
    # amplitude
    x = (π * f₀ * t) .^ 2
    A = (1 .- 2x) .* exp.(-x)
end

function d2t(lyr::Dict{Symbol,Vector{Float64}})
    vₚ = lyr[:vp]
    d = lyr[:dh]
    t = cumsum(2d ./ vₚ)
    pushfirst!(t,0)
end

function d2t(lyr::Vector{Dict{Symbol,Vector{Float64}}})
    stacked = stack(lyr)
    d2t(stacked)
end
        