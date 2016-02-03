##############################################################################
##
## Type
##
##############################################################################

type BansalYaronProblem 

    # consumption process parameters
    μ::Float64 
    νD::Float64
    κμ::Float64 
    κσ::Float64 
    νμ::Float64 
    νσ::Float64 

    # utility parameters
    ρ::Float64  
    γ::Float64 
    ψ::Float64
    θ::Float64

    # algorithm parameters
    invdμ::Float64
    invdσ::Float64 

    # grid    
    μs::Vector{Float64}
    σs::Vector{Float64}

    # initial value
    V::Vector{Float64}
end

function BansalYaronProblem(;μ = 0.0015, νD = 0.0078, κμ = 0.0212, κσ = 0.0131, νμ = 0.0003432, νσ = 0.0378,  ρ = 0.002, γ = 7.5, ψ = 1.5, μn = 50, σn = 50)

    θ = (1-γ)/(1-1/ψ)

    # create drift grid +-3 sd of stationary distribution
    μmin = μ - 6 * νμ / sqrt(2*κμ)
    μmax = μ + 6 * νμ / sqrt(2*κμ)
    μs = collect(linspace(μmin, μmax, μn))
    invdμ = (μn - 1)/(μmax - μmin)

    # create volatility +-3 sd otherwise negative volatility
    σmin = 1.0 - 3 * νσ / sqrt(2*κσ)
    σmax = 1.0 + 3 * νσ / sqrt(2*κσ)
    if σmin < 0
        σmin = 1e-3
        σmax = 2.0
    end
    σs = collect(linspace(σmin, σmax, σn))
    invdσ = (σn - 1)/(σmax - σmin)

    # initial value
    V = Array(Float64, μn * σn)
    fill!(V, (-1/(θ*ρ) * (μ * (1-γ) - 0.5 * (1-γ) * γ * νD^2 * 1.0) + 1.0)^(-1/(1-1/θ)))

    BansalYaronProblem(μ, νD, κμ , κσ , νμ , νσ, ρ, γ, ψ, θ, invdμ, invdσ, μs, σs, V)
end



solve(byp::BansalYaronProblem; method = :trust_region) = solve(Val{method}, byp)

##############################################################################
##
## Plot
##
##############################################################################

function plot(byp::BansalYaronProblem, solution, symbol::Symbol)
    pd = 1/byp.θ * log(max(solution, 0.0)) - log(byp.ρ)
    m = repeat(byp.μs, outer = [length(byp.σs)])
    s2 = repeat(byp.νD^2 * byp.σs, inner = [length(byp.μs)])
    df = DataFrame(wealthconsumption = pd, drift = m, volatility = s2)
    if symbol == :s2
        df = df[df[:, :volatility] .< 0.00010, :]
        condition = find((df[:drift] .== byp.μs[1]) | (df[:drift] .== byp.μs[50]) | (df[:drift] .== byp.μs[end]))
        plot(df[condition, :], x = "volatility", y = "wealthconsumption", color = "drift", Geom.line)
    elseif symbol == :m
        condition = find((df[:volatility] .== byp.νD^2 * byp.σs[1]) | (df[:volatility] .== byp.νD^2 * byp.σs[50]) | (df[:volatility] .== byp.νD^2 * byp.σs[end]))
        plot(df[condition, :], x = "drift", y = "wealthconsumption", color = "volatility", Geom.line)
    end
end

function plot_ll(byp::BansalYaronProblem, symbol::Symbol)
    A1 = (1 - 1 / byp.ψ) / (1 - 0.997 * exp(-byp.κμ))
    A2 = 0.5 * byp.θ * ((1 - 1 / byp.ψ)^2 + (A1 * 0.997 * byp.νμ / byp.νD)^2) / (1 - 0.997 * exp(-byp.κσ))
    m = repeat(byp.μs, outer = [length(byp.σs)])
    s2 = repeat(byp.νD^2 * byp.σs, inner = [length(byp.μs)])
    pd = A1 .* m + A2 .* s2
    df = DataFrame(wealthconsumption = pd, drift = m, volatility = s2)
    if symbol == :s2
        df = df[df[:, :volatility] .< 0.00010, :]
        condition = find((df[:drift] .== byp.μs[1]) | (df[:drift] .== byp.μs[50]) | (df[:drift] .== byp.μs[end]))
        plot(df[condition, :], x = "volatility", y = "wealthconsumption", color = "drift", Geom.line)
    elseif symbol == :m
        condition = find((df[:volatility] .== byp.νD^2 * byp.σs[1]) | (df[:volatility] .== byp.νD^2 * byp.σs[50]) | (df[:volatility] .== byp.νD^2 * byp.σs[end]))
        plot(df[condition, :], x = "drift", y = "wealthconsumption", color = "volatility", Geom.line)
    end
end

