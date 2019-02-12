@recipe function f(p::EmpiricalPattern; variables=[:rho, :vp])
    @assert variables ⊆ names(p) "Variable name(s) not found"
    depth = cumsum(first(reals(p, :dh)))
    l = labels(p)
    
    layout --> (1, length(variables) + 1)
    legend --> false
    yflip --> true
    
    @series begin
        subplot := 1
        seriestype --> :bar
        linewidth --> 0
        ylims --> (0.5, length(l) + 0.5)
        bar_width --> 1.
        fillcolor --> l
        orientation --> :horizontal
        
        l
    end
    
    for (i, var) in enumerate(variables)
        realizations = reals(p, var)
        
        ylims --> (0, depth[end])
        for r in realizations
            @series begin
                subplot := i + 1
                linecolor --> :grey
                r, depth
            end
        end
        
        @series begin
            subplot := i + 1
            linecolor --> :red
            title --> string(var)
            mean(realizations), depth
        end
    end
end