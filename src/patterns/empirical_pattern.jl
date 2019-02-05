struct EmpiricalPattern <: AbstractPattern
    labels::Vector{Int64}
    realizations::Vector{Dict{Symbol, Vector{Float64}}}
    δh::Float64
end

# return realization count
realcount(pattern::EmpiricalPattern) = length(pattern.realizations)

# return realizations
reals(pattern::EmpiricalPattern) = pattern.realizations
reals(pattern::EmpiricalPattern, name::Symbol) = getindex.(reals(p), name)

# sample from realizations
sample(pattern::EmpiricalPattern) = sample(reals(pattern))
sample(pattern::EmpiricalPattern, n::Int64) = sample(reals(patterns), n)

# return base pattern labels
labels(pattern::EmpiricalPattern) = pattern.labels

# show methods
function show(io::IO, p::EmpiricalPattern)
    count = realcount(p)
    print(io, "EmpiricalPattern($count)")
end