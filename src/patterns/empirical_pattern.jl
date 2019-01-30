struct EmpiricalPattern <: AbstractPattern
    labels::Vector{Real}
    properties::Vector{Dict{Symbol, Vector}}
    δh::Float64
end

realcount(pattern::EmpiricalPattern) = length(pattern.properties)

function getreals(pattern::EmpiricalPattern, name::Symbol)
    # extract all realizations of properties
    reals = [vals[name] for vals in pattern.properties]
end

function sample(pattern::EmpiricalPattern) 
    sample(pattern.properties)
    
end

sample(pattern::EmpiricalPattern, n::Int64) = [sample(pattern) for _=1:n]

labels(pattern::EmpiricalPattern) = pattern.labels

properties(pattern::EmpiricalPattern) = pattern.properties

function show(io::IO, p::EmpiricalPattern)
    count = realcount(p)
    print(io, "EmpiricalPattern($count)")
end