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

    σ = sqrt(byp.νμ^2 / (2 * byp.κμ))
    μmin = quantile(Normal(byp.μ, σ), 0.01)
    μmax = quantile(Normal(byp.μ, σ), 0.99)

    α =  2 * byp.κσ / byp.νσ^2
    σmin = quantile(Gamma(α, α), 0.01)
    σmax = quantile(Gamma(α, α), 0.99)
    df = df[(df[:volatility] .> byp.νD^2 * σmin) & (df[:volatility] .< byp.νD^2 * σmax), :]
    df = df[(df[:drift] .> μmin) & (df[:drift] .< μmax), :]
    if symbol == :s2
        rank = denserank(df[:drift])
        condition = (rank .== 1) | (rank .== maximum(rank)) | (rank .== div(maximum(rank), 2))
        convert(Vector{Bool}, map(x -> isapprox(x, minimum(df[:drift])) | isapprox(x, median(df[:drift])) | isapprox(x, maximum(df[:drift])), df[:drift]))
        plot(df[condition, :], x = "volatility", y = "wealthconsumption", color = "drift", Geom.line)
    elseif symbol == :m
        rank = denserank(df[:volatility])
        condition = (rank .== 1) | (rank .== maximum(rank)) | (rank .== div(maximum(rank), 2))
        plot(df[condition, :], x = "drift", y = "wealthconsumption", color = "volatility", Geom.line)
    end
end

