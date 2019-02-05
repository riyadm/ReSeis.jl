struct PatternSet
    patterns::Vector{EmpiricalPattern}
    distributions::Vector{Distribution}
end

function patternset(data::DataFrame, distributions::Vector, w::Int64, s::Int64; 
        flagkey=:facies, propkeys=[:rho, :vp])
    
    @assert all(.!isnan.(convert(Matrix, data[propkeys]))) "NaN values present in properties."
    @assert flagkey in names(data) "Flag not found."
    @assert maximum(data[flagkey]) == length(distributions) "Number of distributions should match the number of classes."
    
    # depth sampling rate
    δh = data[:depth][2]-data[:depth][1]
    
    # vector of unique sequences
    sequences = Dict()
    
    for i in 1:s:size(data, 1)-w
        labels = view(data[flagkey], i:i+w-1)
        seq, counts = seqcounts(labels)
        
        props = Dict(key => view(data[key], i:i+w-1) for key in propkeys)
        props[:thickness] = δh * counts
        props[:dh] = δh * ones(w)
        
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
        patternset(data, distributions, w, s; kwargs...)
end

patterns(pset::PatternSet) = pset.patterns

function sample(pset::PatternSet; weighted=false)
    sample(pset.patterns) # TODO: weighted by frequency
end

sample(pset::PatternSet, n::Int64) = sample(pset.patterns, n)

function sampleresponse(pset::PatternSet, n=3; cut=2, kwargs...)
    @assert 0 <= cut <= n "Pattern to cut doesn't exist."
    @assert in(:dt, keys(kwargs)) "Keyword argument dt missing."
    
    dt = kwargs[:dt]
    
    # sample n patterns from the set
    patterns = sample(pset, n)
    
    # sample property realizations for each pattern
    properties = sample.(patterns)
    
    @show mean([mean(p[:vp]) for p in properties])
    @show mean(properties[cut][:vp])
    # compute synthetic seismogram
    s = synthetic(properties; kwargs...)
    
    if cut > 0
        t = d2t(properties)
     
        # indices to cut
        len = length(first(values(first(properties)))) # TODO: Add length to pattern object or find a better solution to this
        start = (cut-1)*len
        stop = cut*len + 1 
        t_s = dt*(0:(length(s)-1))
        idx = (t_s .>= t[start]) .& (t_s .<= t[stop])
        response = s[idx]
    else
        s
    end    
end

function seqcounts(labels::Vector)
    # find sequence boundaries
    indices = findall(diff(labels) .!= 0)
    
    # get sequence of unique labels
    sequence = labels[indices]
    
    # get label counts
    push!(sequence, labels[end])
    pushfirst!(indices, 0)
    push!(indices, length(labels))
    counts = diff(indices)
    
    sequence, Int64.(counts)
end

function Base.show(io::IO, p::PatternSet)
    count = length(p.patterns)
    numclasses = length(p.distributions)
    print(io, "PatternSet($count, $numclasses)")
end