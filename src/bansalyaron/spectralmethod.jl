
function F!{T}(byp::BansalYaronProblem, μs::Vector{Float64}, σs::Vector{Float64}, bs::Matrix{Float64}, bs∂μ::Matrix{Float64}, bs∂σ::Matrix{Float64}, bs∂2μ::Matrix{Float64}, bs∂2σ::Matrix{Float64}, y::Vector{T}, ydot::Vector{T})
    μn = length(μs)
    σn = length(σs)
    IJ = μn * σn
    ij = 0
    @inbounds for σi in 1:σn, μi in 1:μn
        ij += 1
        V = zero(T)
        @simd for kl in 1:IJ
            V += bs[kl, ij] * y[kl]
        end
        V∂μ = zero(T)
        @simd for kl in 1:IJ
            V∂μ += bs∂μ[kl, ij] * y[kl]
        end
        V∂σ = zero(T)
        @simd for kl in 1:IJ
            V∂σ += bs∂σ[kl, ij] * y[kl]
        end
        V∂2μ = zero(T)
        @simd for kl in 1:IJ
            V∂2μ += bs∂2μ[kl, ij] * y[kl]
        end
        V∂2σ = zero(T)
        @simd for kl in 1:IJ
            V∂2σ += bs∂2σ[kl, ij] * y[kl] 
        end
        if (μi == 1 | μi == μn)
            ydot[ij] = V∂μ
        elseif (σi == 1 | σi == σn)
            ydot[ij] = V∂σ
        else
            ydot[ij] = (byp.ρ * byp.θ * max(zero(T), V^(1-1/byp.θ))
                           + (- byp.ρ * byp.θ + (1-byp.γ) * μs[μi] - 0.5 * (1-byp.γ) * byp.γ * byp.νD^2 * σs[σi]) * V
                           + (byp.κμ * (byp.μ - μs[μi])) * V∂μ
                           + (byp.κσ * (1.0 - σs[σi])) * V∂σ
                           + (0.5 * byp.νμ^2 * σs[σi]) * V∂2μ
                           + (0.5 * byp.νσ^2 * σs[σi]) * V∂2σ
                           )
        end
    end
    return ydot
end

function solve(::Type{Val{:spectral}}, byp::BansalYaronProblem, μs, σs,  bs::Matrix{Float64}, bs∂μ::Matrix{Float64}, bs∂σ::Matrix{Float64}, bs∂2μ::Matrix{Float64}, bs∂2σ::Matrix{Float64}, oldV; kwargs...)
    out = nlsolve((y, ydot) -> F!(byp, μs, σs,  bs, bs∂μ, bs∂σ, bs∂2μ, bs∂2σ, y, ydot), oldV; method = :trust_region, show_trace = true, kwargs...)
    V = bs' * out.zero
    return μs, σs, V
end



function solve(::Type{Val{:spectral}}, byp::BansalYaronProblem; iterations = 50, kwargs...)
    μmin = byp.μ - 6 * byp.νμ / sqrt(2*byp.κμ)
    μmax = byp.μ + 6 * byp.νμ / sqrt(2*byp.κμ)
    σmin = 0.0
    σmax = 3.0
    basisμ = Basis(Cheb, 25, μmin, μmax)
    basisσ = Basis(Cheb, 25, σmin, σmax)
    basis = Basis(basisμ, basisσ)
    S, (μs, σs) = nodes(basis)
    bs = BasisStructure(basis, Expanded(), S, [0 0]').vals[1]'
    bs∂μ = BasisStructure(basis, Expanded(), S, [1 0]').vals[1]'
    bs∂σ = BasisStructure(basis, Expanded(), S, [0 1]').vals[1]'
    bs∂2μ = BasisStructure(basis, Expanded(), S, [2 0]').vals[1]'
    bs∂2σ = BasisStructure(basis, Expanded(), S, [0 2]').vals[1]'
    oldV = bs' \ fill(byp.V, size(bs, 1))
    solve(Val{:spectral}, byp, μs, σs, bs, bs∂μ, bs∂σ, bs∂2μ, bs∂2σ, oldV; iterations = iterations, kwargs...)
end

