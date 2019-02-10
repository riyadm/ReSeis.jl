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

(sim::Kennet)(lyr::Dict{Symbol,Vector{Float64}}) = simkennet(lyr, sim.wavelet, sim.dt, sim.mopt, sim.fs, sim.varnames)
(sim::Kennet)(p::AbstractPattern) = sim(sample(p))
(sim::Kennet)(ps::Vector{<:AbstractPattern}, cut::Int=0) = sim(sample.(ps), cut)

function (sim::Kennet)(lyr::Vector{Dict{Symbol,Vector{Float64}}}, cut::Int=0) 
    s = sim(stack(lyr))

    if cut > 0
        t = d2t(lyr)
        len = length(first(values(first(lyr))))
        start = (cut - 1) * len
        stop = cut * len + 1 
        t_s = sim.dt * (0:(length(s) - 1))
        range = (t_s .>= t[start]) .& (t_s .<= t[stop])
        s = s[range]
    end
    
    s
end


function (sim::Kennet)(pset::PatternSet; mask::Vector{Int}=[0, 0, 0], cut::Int=2)
    n = length(mask)
    0 <= cut <= n || error("Attempt to cut $n-element stack at index $cut")
    
    # sample n patterns from the set
    idx = map(x -> x == 0 ? rand(1:length(pset)) : x, mask)
    patterns = pset.patterns[idx]
    
    # compute synthetic seismogram
    s = sim(patterns, cut)
    
    s
end

function simkennet(lyr::Dict{Symbol,Vector{Float64}}, wavelet::Vector{Float64}, dt::Float64=1e-3, mopt::Int64=2, fs::Int64=0, 
        varnames::Vector{Symbol}=[:rho, :vp, :dh])
    
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

function Base.show(io::IO, sim::Kennet)
    print(io, "Kennet($sim.varnames)")
end

function Base.show(io::IO, ::MIME"text/plain", sim::Kennet)
    println(io, "  Sampling rate Δt:  ", sim.dt, " s")
    println(io, "  Variables:         ", sim.varnames)
    println(io, "  Multiples order:   ", sim.mopt) 
    println(io, "  Wavelet length:    ", length(sim.wavelet), " points")
end
