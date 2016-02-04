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

function BansalYaronFDMethod(byp::BansalYaronProblem; μn = 50, σn = 50)
    # create drift grid +-3 sd of stationary distribution
    μmin = byp.μ - 6 * byp.νμ / sqrt(2*byp.κμ)
    μmax = byp.μ + 6 * byp.νμ / sqrt(2*byp.κμ)
    μs = collect(linspace(μmin, μmax, μn))
    invdμ = (μn - 1)/(μmax - μmin)
  
    σmin = 0.0
    σmax = 3.0
    σs = collect(linspace(σmin, σmax, σn))
    invdσ = (σn - 1)/(σmax - σmin)

    BansalYaronFDMethod(byp, invdμ, invdσ, μs, σs)
end

function solve(x::Union{Type{Val{:ode}}, Type{Val{:nl}}}, byp::BansalYaronProblem; μn = 50, σn = 50,  kwargs...)
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
        if σi == 1
            V∂σ = - 0.5 * (3 * fy[μi, σi] - 4 * fy[μi, σi+1] + fy[μi, σi+2]) * byfd.invdσ
        elseif σi == σn
            V∂σ = 0.5 * (3 * fy[μi, σi] - 4 * fy[μi, σi-1] + fy[μi, σi-2]) * byfd.invdσ
        else
            V∂2σ = (fy[μi, σi+1] + fy[μi, σi-1] - 2.0 * fy[μi, σi]) * byfd.invdσ^2
            V∂σ = 0.5 * (fy[μi, σi+1] - fy[μi, σi-1]) * byfd.invdσ
        end
        V∂μ = zero(T)
        V∂2μ = zero(T)
        if μi == 1
            V∂μ = - 0.5 * (3 * fy[μi, σi] - 4 * fy[μi+1, σi] + fy[μi+2, σi]) * byfd.invdμ
        elseif μi == μn
            V∂μ = 0.5 * (3 * fy[μi, σi] - 4 * fy[μi-1, σi] + fy[μi-2, σi]) * byfd.invdμ
        else
            V∂2μ = (fy[μi+1, σi] + fy[μi-1, σi] - 2.0 * fy[μi, σi]) * byfd.invdμ^2
            V∂μ = 0.5 * (fy[μi+1, σi] - fy[μi-1, σi]) * byfd.invdμ
        end
        ydot[ij] = (byp.ρ * byp.θ * max(V, 0.0)^(1-1/byp.θ)
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