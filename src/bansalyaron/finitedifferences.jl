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
    σ = sqrt(byp.νσ^2 / (2 * byp.κσ))
    σmin = max(0.01, quantile(Normal(1.0, σ), 0.001))
    σmax = quantile(Normal(1.0, σ), 0.999)
    σs = collect(linspace(σmin, σmax, σn))
    invdσ = 1 / (σs[2] - σs[1])

    # create drift grid
    σ = sqrt(byp.νμ^2 / (2 * byp.κμ))
    μmin = quantile(Normal(byp.μ, σ), 0.001)
    μmax = quantile(Normal(byp.μ, σ), 0.999)
    μs = collect(linspace(μmin, μmax, μn))
    invdμ = 1 / (μs[2] - μs[1])

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
function F!{T}(byfd::BansalYaronFDMethod, y::Vector{T}, ydot::Vector{T})
    byp = byfd.byp ; μs = byfd.μs ; σs = byfd.σs
    μn = length(μs)
    σn = length(σs)
    fy = reshape(y, μn, σn)
    ij = 0
    for σi in 1:σn, μi in 1:μn
        ij += 1
        G = fy[μi, σi]
        G∂σ = zero(T)
        G∂2σ = zero(T)
        if (byp.κσ * (1.0 - σs[σi])) >= 0
            G∂σ = (fy[μi, σi+1] - fy[μi, σi]) * byfd.invdσ
        else
            G∂σ = (fy[μi, σi] - fy[μi, σi-1]) * byfd.invdσ
        end
        if σi == 1
            G∂2σ = (fy[μi, σi+1] - fy[μi, σi]) * byfd.invdσ^2
        elseif σi == σn
            G∂2σ = (fy[μi, σi-1] - fy[μi, σi]) * byfd.invdσ^2
        else
            G∂2σ = (fy[μi, σi+1] + fy[μi, σi-1] - 2 * fy[μi, σi]) * byfd.invdσ^2
        end
        G∂μ = zero(T)
        G∂2μ = zero(T)
        if (byp.κμ * (byp.μ - μs[μi])) >= 0
            G∂μ = (fy[μi+1, σi] - fy[μi, σi]) * byfd.invdμ
        else
            G∂μ = (fy[μi, σi] - fy[μi-1, σi]) * byfd.invdμ
        end
        if μi == 1
            G∂2μ = (fy[μi+1, σi] - fy[μi, σi]) * byfd.invdμ^2
        elseif μi == μn
            G∂2μ = (fy[μi-1, σi] - fy[μi, σi]) * byfd.invdμ^2
        else
            G∂2μ = (fy[μi+1, σi] + fy[μi-1, σi] - 2 * fy[μi, σi]) * byfd.invdμ^2
        end
        # note G / byp.θ for method of line
        ydot[ij] = G / byp.θ * (byp.θ / G - byp.θ * byp.ρ + (1 - byp.γ) * (μs[μi] - 0.5 * byp.γ * byp.νD^2 * σs[σi])
                  + byp.θ  * (byp.κμ * (byp.μ - μs[μi])) * G∂μ / G 
                  + byp.θ  * (byp.κσ * (1.0 - σs[σi])) * G∂σ / G 
                  + 0.5 * byp.θ * byp.νμ^2 * σs[σi] * (G∂2μ / G + (byp.θ - 1) * G∂μ^2 / G^2)
                  + 0.5 * byp.θ * byp.νσ^2 * (G∂2σ / G + (byp.θ - 1) * G∂σ^2 / G^2)
                  )
    end
    return ydot
end

function solve(::Type{Val{:ode}}, byfd::BansalYaronFDMethod; maxiter = 100, tol = 1e-5, kwargs...)
    oldG = fill(byfd.byp.G, length(byfd.μs) * length(byfd.σs))
    newG = oldG
    output = deepcopy(newG)
    distance = Inf
    i = 0
    ydot = deepcopy(oldG)
    while i < maxiter
        i += 1
        xout, yout = ODE.ode45((t, y) -> (ydot = deepcopy(y); F!(byfd, y, ydot)), oldG, [0.0 ; 10.0]; minstep = 1e-3, kwargs...)
        newG = yout[end]
        distance =  maxabs(F!(byfd, newG, output))
        @show i, distance
        if distance < tol
            break
        else
            oldG = newG
        end
    end
    return newG
end


function solve(::Type{Val{:nl}}, byfd::BansalYaronFDMethod; maxiter = 100, tol = 1e-5, kwargs...)
    oldG = fill(byfd.byp.G, length(byfd.μs) * length(byfd.σs))
    out = nlsolve((y, ydot) -> F!(byfd, y, ydot), oldG, method = :trust_region, show_trace = true, autodiff = true, iterations = maxiter, ftol = tol, kwargs...)
    return out.zero
end