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
        labels = data[flagkey][i:i+w-1]
        seq, counts = seqcounts(labels)
        
        props = Dict(key => data[key][i:i+w-1] for key in propkeys)
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

function sample(pset::PatternSet; weighted=false)
    sample(pset.patterns) # TODO: weighted by frequency
end

sample(pset::PatternSet, n::Int64) = [sample(pset) for _=1:n]

function sampleresponse(pset::PatternSet, n=3; cut=2, kwargs...)
    @assert cut <= n "Pattern to cut doesn't exist."
    
    # sample n patterns from the set
    patterns = sample(pset, n)
    
    # sample property realizations for each pattern
    properties = sample.(patterns)

    s = synthetic(properties; kwargs...)
    t = d2t(properties)
    
    # indices to cut
    len = length(first(values(first(properties)))) # TODO: Add length to pattern object or find a better solution to this
    start = (cut-1)*len
    stop = cut*len + 1
    dt = kwargs[:dt]
    t_s = collect(dt*(0:(length(s)-1)))
    @show length(t_s), length(s)
    response = s[(t_s .> t[start]) .& (t_s .< t[stop])]
    
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