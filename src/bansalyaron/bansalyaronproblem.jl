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

    # initial value
    V::Float64
end

function BansalYaronProblem(;μ = 0.0015, νD = 0.0078, κμ = 0.0212, κσ = 0.0131, νμ = 0.0003432, νσ = 0.0378,  ρ = 0.002, γ = 7.5, ψ = 1.5)

    θ = (1-γ)/(1-1/ψ)

    # initial value
    V = (-1/(θ*ρ) * (μ * (1-γ) - 0.5 * (1-γ) * γ * νD^2 * 1.0) + 1.0)^(-1/(1-1/θ))

    BansalYaronProblem(μ, νD, κμ , κσ , νμ , νσ, ρ, γ, ψ, θ, V)
end



solve(byp::BansalYaronProblem; method = :trust_region) = solve(Val{method}, byp)

##############################################################################
##
## Plot
##
##############################################################################

function plot(byp::BansalYaronProblem, solution, symbol::Symbol)
    μs, σs, V = solution
    pd = 1/byp.θ * log(V) - log(byp.ρ)
    m = repeat(μs, outer = [length(σs)])
    s2 = repeat(byp.νD^2 * σs, inner = [length(μs)])
    df = DataFrame(wealthconsumption = pd, drift = m, volatility = s2)
    if symbol == :s2
        df = df[df[:, :volatility] .< 0.00010, :]
        condition = find((df[:drift] .== μs[1]) | (df[:drift] .== μs[div(length(μs), 2)]) | (df[:drift] .== μs[end]))
        plot(df[condition, :], x = "volatility", y = "wealthconsumption", color = "drift", Geom.line)
    elseif symbol == :m
        condition = find((df[:volatility] .== byp.νD^2 * σs[1]) | (df[:volatility] .== byp.νD^2 * σs[div(length(σs), 2)]) | (df[:volatility] .== byp.νD^2 * σs[end]))
        plot(df[condition, :], x = "drift", y = "wealthconsumption", color = "volatility", Geom.line)
    end
end

