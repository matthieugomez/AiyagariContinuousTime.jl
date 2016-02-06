##############################################################################
##
## Finite Difference
##
##############################################################################

type BansalYaronFDMethod
    byp::BansalYaronProblem
    
    # algorithm parameters
    invdμ::Float64
    invdσ::Float64 

    # grid    
    μs::Vector{Float64}
    σs::Vector{Float64}
end

function BansalYaronFDMethod(byp::BansalYaronProblem; μn = 30, σn = 50)

    # create volatility grid. square root process has stationary Gamma distribution
    σmin = 0.0
    α =  2 * byp.κσ / byp.νσ^2
    σmax = quantile(Gamma(α, α), 0.999)
    σs = collect(linspace(σmin, σmax, σn))
    invdσ = (σn - 1)/(σmax - σmin)

    # create drift grid
    σ = sqrt(byp.νμ^2 / (2 * byp.κμ))
    μmin = quantile(Normal(byp.μ, σ), 0.001)
    μmax = quantile(Normal(byp.μ, σ), 0.999)
    μs = collect(linspace(μmin, μmax, μn))
    invdμ = (μn - 1)/(μmax - μmin)

    BansalYaronFDMethod(byp, invdμ, invdσ, μs, σs)
end

function solve(x::Union{Type{Val{:ode}}, Type{Val{:nl}}}, byp::BansalYaronProblem; μn = 30, σn = 50,  kwargs...)
    byfd = BansalYaronFDMethod(byp, μn = μn, σn = σn)
    solution = solve(x, byfd, kwargs...)
    return byfd.μs, byfd.σs, solution
end

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

function F!{T}(byfd::BansalYaronFDMethod, y::Vector{T}, ydot::Vector{T})
    byp = byfd.byp ; μs = byfd.μs ; σs = byfd.σs
    μn = length(μs)
    σn = length(σs)
    fy = Grid(y, μn, σn)
    ij = 0
    for σi in 1:σn, μi in 1:μn
        ij += 1
        V = fy[μi, σi]
        V∂σ = zero(T)
        V∂2σ = zero(T)
        if (byp.κσ * (1.0 - σs[σi])) >= 0
            V∂σ = (fy[μi, σi+1] - fy[μi, σi]) * byfd.invdσ
        else
            V∂σ = (fy[μi, σi] - fy[μi, σi-1]) * byfd.invdσ
        end
        if σi == 1
            # does not matter since volatility == 0
            V∂2σ = zero(T)
        elseif σi == σn
            V∂2σ = (fy[μi, σi-1] - fy[μi, σi]) * byfd.invdσ^2
        else
            V∂2σ = (fy[μi, σi+1] + fy[μi, σi-1] - 2 * fy[μi, σi]) * byfd.invdσ^2
        end
        V∂μ = zero(T)
        V∂2μ = zero(T)
        if (byp.κμ * (byp.μ - μs[μi])) >= 0
            V∂μ = (fy[μi+1, σi] - fy[μi, σi]) * byfd.invdμ
        else
            V∂μ = (fy[μi, σi] - fy[μi-1, σi]) * byfd.invdμ
        end
        if μi == 1
            V∂2μ = (fy[μi+1, σi] - fy[μi, σi]) * byfd.invdμ^2
        elseif μi == μn
            V∂2μ = (fy[μi-1, σi] - fy[μi, σi]) * byfd.invdμ^2
        else
            V∂2μ = (fy[μi+1, σi] + fy[μi-1, σi] - 2 * fy[μi, σi]) * byfd.invdμ^2
        end
        ydot[ij] = (byp.ρ * byp.θ * max(V, zero(T))^(1-1/byp.θ)
                    + (- byp.ρ * byp.θ + (1-byp.γ) * μs[μi] - 0.5 * (1-byp.γ) * byp.γ * byp.νD^2 * σs[σi]) * V
                    + (byp.κμ * (byp.μ - μs[μi])) * V∂μ 
                    + (byp.κσ * (1.0 - σs[σi])) * V∂σ 
                    + (0.5 * byp.νμ^2 * σs[σi]) * V∂2μ 
                    + (0.5 * byp.νσ^2 * σs[σi]) * V∂2σ 
                    )
    end
    return ydot
end

function solve(::Type{Val{:ode}}, byfd::BansalYaronFDMethod; iterations = 1000, kwargs...)
    oldV = fill(byfd.byp.V, length(byfd.μs) * length(byfd.σs))
    newV = oldV
    distance = Inf
    i = 0
    ydot = deepcopy(oldV)
    while i < iterations
        i += 1
        xout, yout = ODE.ode45((t, y) -> F!(byfd, y, ydot), oldV, [0.0 ; 100.0]; kwargs...)
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


function solve(::Type{Val{:nl}}, byfd::BansalYaronFDMethod; kwargs...)
    oldV = fill(byfd.byp.V, length(byfd.μs) * length(byfd.σs))
    out = nlsolve((y, ydot) -> F!(byfd, y, ydot), oldV, method = :trust_region, show_trace = true, autodiff = true, kwargs...)
    return out.zero
end