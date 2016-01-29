##############################################################################
##
## Solve
##
##############################################################################

function shift(x::Vector, i, s)
    newi = i + s
    if (newi <= length(x)) & (newi >= 1)
        return x[newi]
    else
        return x[i]
    end
end

function F!(byp::BansalYaronProblem, y::Vector, ydot::Vector)
    ij = zero(Int)
    μn = length(byp.μs)
    σn = length(byp.σs)
    @inbounds for σi in 1:σn, μi in 1:μn
        ij += 1
        ∂μ = byp.κμ * (byp.μ - byp.μs[μi]) * byp.invdμ
        ∂2μ = 0.5 * byp.νμ^2 * byp.σs[σi] * byp.invdμ^2
        ∂σ = byp.κσ * (1.0 - byp.σs[σi]) * byp.invdσ
        ∂2σ = 0.5 * byp.νσ^2 * byp.σs[σi] * byp.invdσ^2
        current = 0.0
        current += ∂2μ * (shift(y, ij, 1) + shift(y, ij, -1) - 2 * shift(y, ij, 0)) 
        if ∂μ > 0 
            current += ∂μ * (shift(y, ij, 1) - shift(y, ij, 0))
        else 
            current += ∂μ * (shift(y, ij, 0) - shift(y, ij, -1))
        end
        current += ∂2σ * (shift(y, ij, μn) + shift(y, ij, -μn) - 2 * shift(y, ij, 0))
        if ∂σ > 0 
            current += ∂σ * (shift(y, ij, μn) - shift(y, ij, 0))
        else 
            current += ∂σ * (shift(y, ij, 0) - shift(y, ij, - μn))
        end 
        current += byp.ρ * byp.θ * max(y[ij], 0.0)^(1-1/byp.θ) + (- byp.ρ * byp.θ + (1-byp.γ) * byp.μs[μi] - 0.5 * (1-byp.γ) * byp.γ * byp.νD^2 * byp.σs[σi])*y[ij] 
        ydot[ij] = current
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
