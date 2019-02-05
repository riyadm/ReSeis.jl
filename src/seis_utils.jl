function simkennet(lyr::Dict{Symbol,Vector{Float64}}; wavelet::Vector=[], dt::Float64=1e-3, 
        varnames::Vector{Symbol}=[:rho, :vp, :dh], mopt::Int64=2, fs::Int64=0)
    
    @assert all(haskey(lyr, key) for key in varnames) "Property name not found."
    @assert wavelet != nothing "No wavelet provided."
    
    
    ρ = lyr[varnames[1]]
    vₚ = lyr[varnames[2]]
    d = lyr[varnames[3]]
    
    nlr = length(ρ)
    n = length(wavelet)
    
    om = (2π/(n*dt)) .* (0:(n/2 - 1.)) |> collect
    freq = om ./ (2π)
    
    p₀ = ifft(wavelet)
    w₀ = (1 ./ (ρ[1]*vₚ[1])) .* p₀
    rdhat = zeros(size(om))
    # tdhat = ones(size(om))
    
    deno = @. ρ[2:end]*vₚ[2:end]+ ρ[1:end-1]*vₚ[1:end-1]
    rd = .- diff(ρ.*vₚ) ./ deno
    td = @. 2*sqrt(ρ[2:end]*vₚ[2:end]*ρ[1:end-1]*vₚ[1:end-1])/deno
    
    if fs == 1
        pushfirst!(rd, -1)
        pushfirst!(td, 1)
    else
        pushfirst!(rd, 0)
        pushfirst!(td, 1)
    end
    
    ru = .- rd
    # tu = td
    
    for j=nlr:-1:1
        ed = exp.((d[j]/vₚ[j])*im*om)
        
        if mopt == 0
            reverb = ones(size(om))
        elseif mopt == 1
            reverb = 1 .+ ru[j] .*ed .*rdhat .*ed
        else
            reverb = 1 ./(1 .- ru[j] .*ed .*rdhat .*ed)
        end
        
        rdhat = rd[j] .+ td[j] .*ed .*rdhat .*ed .*reverb .*td[j]
        #tdhat = tdhat .*ed .*reverb .*td[j]
        
        if mopt == 1
            dx=findall(real(rdhat) .> 1); rdhat[dx] = 1 + im*imag(rdhat[dx]);
            dx=findall(imag(rdhat) .> 1); rdhat[dx] = real(rdhat[dx]) + im;
            # dx=findall(real(tdhat) .> 1); tdhat[dx] = 1 + im*imag(tdhat[dx]);
            # dx=findall(imag(tdhat) .> 1); tdhat[dx] = real(tdhat[dx]) + im;
        end
    end
    
    # tf = hcat(freq, rdhat, tdhat)
    # pz = tdhat .* p₀[1:length(om)]
    wz = rdhat .* p₀[1:length(om)]  
    
    # pz = vcat(pz, 0, conj(reverse(pz[2:end], dims=1)))
    wz = vcat(0, wz[2:end], 0, conj(reverse(wz[2:end], dims=1)))
    
    # pz = real(fft(pz))
    wz = real(fft(wz))
end

function synthetic(lyr::Dict{Symbol,Vector{Float64}}; kwargs...)
    s = simkennet(lyr; kwargs...)
end

function synthetic(lyr::Vector{Dict{Symbol,Vector{Float64}}}; kwargs...)
    stacked = stacklayers(lyr)
    s = synthetic(stacked; kwargs...)
end

function stacklayers(lyr::Vector{Dict{Symbol,Vector{Float64}}})
    stacked = eltype(lyr)()
    
    for key in keys(first(lyr))
        properties = getindex.(lyr, key)
        stacked[key] = vcat(properties...)
    end
    stacked
end


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
    stacked = stacklayers(lyr)
    d2t(stacked)
end
        