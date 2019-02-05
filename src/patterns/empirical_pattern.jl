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

# show methods
function show(io::IO, p::EmpiricalPattern)
    count = realcount(p)
    print(io, "EmpiricalPattern($count)")
end