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

patterns(pset::PatternSet) = pset.patterns
Base.size(pset::PatternSet) = length(pset.patterns)
Base.getindex(pset::PatternSet, i) = getindex(pset.patterns, i)
Base.iterate(pset::PatternSet) = iterate(pset.patterns)
Base.lastindex(pset::PatternSet) = lastindex(pset.patterns)


function sample(pset::PatternSet; weighted=false)
    sample(pset.patterns) # TODO: weighted by frequency
end

sample(pset::PatternSet, n::Int64) = sample(pset.patterns, n)


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