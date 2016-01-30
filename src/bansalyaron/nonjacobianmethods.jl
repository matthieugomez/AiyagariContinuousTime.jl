##############################################################################
##
## Solve
##
##############################################################################
type Grid{T}
    x::Vector{T}
    I::Int
    J::Int
end

function Base.getindex(x::Grid, i, j)
    i1 = min(max(i, 1), x.I)
    j1 = min(max(j, 1), x.J)
    return x.x[i1 + x.I * (j1 - 1)]
end


function F!(byp::BansalYaronProblem, y::Vector, ydot::Vector)
    μn = length(byp.μs)
    σn = length(byp.σs)
    fy = Grid(y, μn, σn)
    ij = 0
    for σi in 1:σn, μi in 1:μn
        ij += 1
        ∂μ = byp.κμ * (byp.μ - byp.μs[μi]) * byp.invdμ
        ∂2μ = 0.5 * byp.νμ^2 * byp.σs[σi] * byp.invdμ^2
        ∂σ = byp.κσ * (1.0 - byp.σs[σi]) * byp.invdσ
        ∂2σ = 0.5 * byp.νσ^2 * byp.σs[σi] * byp.invdσ^2
        
        ydot[ij] = (byp.ρ * byp.θ * max(fy[μi, σi], 0.0)^(1-1/byp.θ)
                    + (- byp.ρ * byp.θ + (1-byp.γ) * byp.μs[μi] - 0.5 * (1-byp.γ) * byp.γ * byp.νD^2 * byp.σs[σi]) * fy[μi, σi] 
                    + max(∂μ, 0.0) * (fy[μi+1, σi] - fy[μi, σi]) 
                    + min(∂μ, 0.0) * (fy[μi, σi] - fy[μi-1, σi])
                    + max(∂σ, 0.0) * (fy[μi, σi+1] - fy[μi, σi]) 
                    + min(∂σ, 0.0) * (fy[μi, σi] -  fy[μi, σi-1])
                    + ∂2μ * (fy[μi+1, σi] + fy[μi-1, σi]  - 2.0 * fy[μi, σi]) 
                    + ∂2σ * (fy[μi, σi+1] + fy[μi, σi-1] - 2.0 * fy[μi, σi]) 
                    )
    end
    @show norm(ydot)
    return ydot
end

function F(byp::BansalYaronProblem, y::Vector)
    ydot = deepcopy(y)
    F!(byp, y, ydot)
end

function solve(::Type{Val{:ode23}}, byp::BansalYaronProblem; kwargs...)
    oldV = deepcopy(byp.V)
    for iter in 1:5
       tout, yout = ODE.ode23((t, y) -> F(byp, y), oldV, [0.0;100.0]; kwargs...)
       newV = yout[end]
       distance = chebyshev(vec(oldV), vec(newV))
       if distance < 1e-10
            return newV
        else
            @show iter distance
            oldV, newV = newV, oldV
        end
    end
    return oldV
end
