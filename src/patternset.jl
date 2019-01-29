struct PatternSet
    patterns::Vector{EmpiricalPattern}
    distributions::Vector{Distribution}
end

function patternset(data::DataFrame, distributions::Vector{Distribution}, w::Int64, s::Int64; 
        flagkey=:gmmfacies, propkeys=[:rho, :vp])
    
    @assert all(.!isnan.(convert(Matrix, data[propkeys]))) "NaN values present in properties."
    @assert flagkey in names(data) "Flag not found."
    
    # depth sampling rate
    δh = data[:depth][2]-data[:depth][1]
    
    # vector of unique sequences
    sequences = Dict()
    
    for i in 1:s:size(data, 1)-w
        labels = data[flagkey][i:i+w-1]
        seq, counts = seqcounts(labels)
        
        props = Dict(key => data[key][i:i+w-1] for key in propkeys)
        props[:thickness] = δh * counts
        
        if seq in keys(sequences)
            push!(sequences[seq], props)
        else
            sequences[seq] = [props]
        end
    end
    
    patterns = [EmpiricalPattern(seq, props, δh) for (seq, props) in sequences]
    counts = realcount.(patterns)
    
    # sort patterns in descending order of realization counts
    permute!(patterns, sortperm(counts, rev=true))
        
    PatternSet(patterns, distributions)
end

function patternset(data::DataFrame, gmm::GMM, w::Int64, s::Int64; kwargs...)
        distributions = [MvNormal(means(gmm)[i,:], covars(gmm)[i]) for i=1:gmm.n]
        patternset(data, distributions, w, s, kwargs...)
end

function sample(pset::PatternSet)
    sample(pset.patterns) # weighted by frequency
end

function sample(pset::PatternSet, n::Int64)
    [sample(pset) for _=1:n]
end

function seqcounts(labels::Vector)
    # find sequence boundaries
    indices = findall(labels[1:end-1] .!= labels[2:end])
    sequence = labels[indices]
    push!(sequence, labels[end])
    pushfirst!(indices, 0)
    push!(indices, length(labels))
    counts = diff(indices)
    
    sequence, Int64.(counts)
end

function show(io::IO, p::PatternSet)
    count = length(p.patterns)
    numclasses = length(p.distributions)
    print(io, "PatternSet($count, $numclasses)")
end