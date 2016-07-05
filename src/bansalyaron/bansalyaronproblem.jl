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
    G::Float64
end

function BansalYaronProblem(;μ = 0.018, νD = 0.025, κμ = 0.3, κσ = 0.012, νμ = 0.0114, νσ = 0.189,  ρ = 0.0132, γ = 7.5, ψ = 1.5)
    θ = (1 - γ) / (1 - 1 / ψ)
    G = 1 / (ρ - (1 - 1 / ψ) * (μ - γ / 2 * νD^2))
    BansalYaronProblem(μ, νD, κμ , κσ , νμ , νσ, ρ, γ, ψ, θ, G)
end



solve(byp::BansalYaronProblem; method = :trust_region, kwargs...) = solve(Val{method}, byp ; kwargs...)

##############################################################################
##
## Plot
##
##############################################################################

function plot(byp::BansalYaronProblem, solution, symbol::Symbol)
    μs, σs, G = solution
    m = repeat(μs, outer = [length(σs)])
    s2 = repeat(byp.νD^2 * σs, inner = [length(μs)])
    df = DataFrame(wealthconsumption = log(G), drift = m, volatility = s2)

    
    σ = sqrt(byp.νμ^2 / (2 * byp.κμ))
    μmin = quantile(Normal(byp.μ, σ), 0.01)
    μmax = quantile(Normal(byp.μ, σ), 0.99)

    σ = sqrt(byp.νσ^2 / (2 * byp.κσ))
    σmin = max(0.0, quantile(Normal(1.0, σ), 0.01))
    σmax = quantile(Normal(1.0, σ), 0.99)
    df = df[(df[:volatility] .> byp.νD^2 * σmin^2) & (df[:volatility] .< byp.νD^2 * σmax^2), :]
    df = df[(df[:drift] .> μmin) & (df[:drift] .< μmax), :]

    # using Plotly
    # μ, σ = unique(df[:drift]), sqrt(unique(df[:volatility]))
    # trace = surface(x =μ, y = σ, z = reshape(exp(df[:wealthconsumption]), length(μ), length(σ)))
    # layout = Plotly.Layout(title="Long Run Risk (BKY 2007)", autosize=false,  paper_bgcolor= "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)", scene = attr(xaxis = attr(title = "μ"), yaxis = attr(title = "σ"), zaxis = attr(title = "Wealth / Consumption")))
    # myplot = Plotly.plot(trace, layout)

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

