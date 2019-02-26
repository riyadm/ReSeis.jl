@recipe function f(p::EmpiricalPattern; variables=[:rho, :vp], grad=:rainbow)
    @assert variables ⊆ names(p) "Variable name(s) not found"
    depth = cumsum(first(reals(p, :dh)))
    l = labels(p)
    units = ["kg/m³", "ft/s"]
    layout --> (1, length(variables) + 1)
    legend --> false
    yflip --> true
    size --> (800,400)
    @series begin
        subplot := 1
        seriestype --> :bar
        linewidth --> 0
        ylims --> (0.5, length(l) + 0.5)
        bar_width --> 1.
        seriescolor --> grad
        clims --> (1, 14)
        fill_z --> l
        orientation --> :horizontal
        yticks --> 1:length(l)
        xticks --> 1:14
        ylabel --> "Layer Index"
        xlabel --> "Flag"
        
        l
    end
    
    for (i, var) in enumerate(variables)
        realizations = reals(p, var)
        
        ylims --> (0, depth[end])
        ylabel --> "Depth, ft"
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
            xlabel --> units[i]
            mean(realizations), depth
        end
        
        
        
    end
end