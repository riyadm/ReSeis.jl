struct EmpiricalPattern <: AbstractPattern
    labels::Vector{Int64}
    realizations::Vector{Dict{Symbol, Vector{Float64}}}
    δh::Float64
end

# return realization count
realcount(p::EmpiricalPattern) = length(p.realizations)

# return realizations
reals(p::EmpiricalPattern) = p.realizations
reals(p::EmpiricalPattern, name::Symbol) = getindex.(reals(p), name)

# sample from realizations
sample(p::EmpiricalPattern) = sample(reals(p))
sample(p::EmpiricalPattern, n::Int64) = sample(reals(p), n)

# return base pattern labels
labels(p::EmpiricalPattern) = p.labels

# variable names
Base.names(p::EmpiricalPattern) = collect(keys(first(reals(p))))

# show methods
function Base.show(io::IO, p::EmpiricalPattern)
    print(io, "EmpiricalPattern($(join(labels(p),"-")))")
end

function Base.show(io::IO, ::MIME"text/plain", p::EmpiricalPattern)
    println(io, p)
    println(io, "  Sampling rate Δh:  ", p.δh)
    println(io, "  Variables:         ", join(names(p), ", "))
    println(io, "  Realization count: ", realcount(p)) 
end

