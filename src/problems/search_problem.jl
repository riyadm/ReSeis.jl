struct SearchProblem <: AbstractProblem
    patternset::PatternSet
    trace::Vector{Float64}
    simulator::AbstractSimulator
    solver::AbstractSolver
end

