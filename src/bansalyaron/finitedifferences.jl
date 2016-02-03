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
Base.getindex(g::Grid, i, j) = g.x[i + g.I * (j - 1)]

function F!{T}(byp::BansalYaronProblem, y::Vector{T}, ydot::Vector{T})
    μn = length(byp.μs)
    σn = length(byp.σs)
    fy = Grid(y, μn, σn)
    ij = 0
    for σi in 1:σn, μi in 1:μn
        ij += 1
        V∂σ = zero(T)
        V∂2σ = zero(T)
        if  (σi > 1) & (σi < σn)
            V∂2σ = (fy[μi, σi+1] + fy[μi, σi-1] - 2.0 * fy[μi, σi]) 
            V∂σ = 0.5 * (fy[μi, σi+1] - fy[μi, σi-1]) 
        elseif σi == 1
            V∂σ = - 0.5 * (3 * fy[μi, σi] - 4 * fy[μi, σi+1] + fy[μi, σi+2])
        elseif σi == σn
            V∂σ = 0.5 * (3 * fy[μi, σi] - 4 * fy[μi, σi-1] + fy[μi, σi-2])
        end
        V∂μ = zero(T)
        V∂2μ = zero(T)
        if  (μi > 1) & (μi < μn)
            V∂2μ = (fy[μi+1, σi] + fy[μi-1, σi] - 2.0 * fy[μi, σi]) 
            V∂μ = 0.5 * (fy[μi+1, σi] - fy[μi-1, σi]) 
        elseif μi == 1
            V∂μ = - 0.5 * (3 * fy[μi, σi] - 4 * fy[μi+1, σi] + fy[μi+2, σi])
        elseif μi == μn
            V∂μ = 0.5 * (3 * fy[μi, σi] - 4 * fy[μi-1, σi] + fy[μi-2, σi])
        end
        ydot[ij] = (byp.ρ * byp.θ * max(fy[μi, σi], 0.0)^(1-1/byp.θ)
                    + (- byp.ρ * byp.θ + (1-byp.γ) * byp.μs[μi] - 0.5 * (1-byp.γ) * byp.γ * byp.νD^2 * byp.σs[σi]) * fy[μi, σi] 
                    + (byp.κμ * (byp.μ - byp.μs[μi]) * byp.invdμ) * V∂μ 
                    + (byp.κσ * (1.0 - byp.σs[σi]) * byp.invdσ) * V∂σ 
                    + (0.5 * byp.νμ^2 * byp.σs[σi] * byp.invdμ^2) * V∂2μ 
                    + (0.5 * byp.νσ^2 * byp.σs[σi] * byp.invdσ^2) * V∂2σ 
                    )
    end
    return ydot
end

function F(byp::BansalYaronProblem, y::Vector)
    ydot = deepcopy(y)
    F!(byp, y, ydot)
end

function solve(::Type{Val{:ode}}, byp::BansalYaronProblem; iterations = 1000, kwargs...)
    oldV = deepcopy(byp.V)
    newV = oldV
    distance = Inf
    i = 0
    ydot = deepcopy(oldV)
    while i < iterations
        i += 1
        xout, yout = ODE.ode45((t, y) -> F!(byp, y, ydot), oldV, [0.0 ; 100.0]; kwargs...)
        newV = yout[end]
        distance = chebyshev(newV, oldV)
        @show i, distance
        if distance < 1e-5
            break
        else
            oldV = newV
        end
    end
    return newV
end


function solve(::Type{Val{:nl}}, byp::BansalYaronProblem; kwargs...)
    out = nlsolve((y, ydot) -> F!(byp, y, ydot), byp.V, method = :trust_region, show_trace = true, autodiff = true, kwargs...)
    return out.zero
end