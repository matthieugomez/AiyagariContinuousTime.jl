
function F!(byp::BansalYaronProblem, μs::Vector{Float64}, σs::Vector{Float64}, bs::Matrix{Float64}, bs∂μ::Matrix{Float64}, bs∂σ::Matrix{Float64}, bs∂2μ::Matrix{Float64}, bs∂2σ::Matrix{Float64}, y, ydot)
    μn = length(μs)
    σn = length(σs)
    IJ = μn * σn
    ij = 0
    @inbounds for σi in 1:σn, μi in 1:μn
        ij += 1
        V = zero(eltype(y))
        @simd for kl in 1:IJ
            V += bs[kl, ij] * y[kl]
        end
        V∂μ = zero(eltype(y))
        @simd for kl in 1:IJ
            V∂μ += bs∂μ[kl, ij] * y[kl]
        end
        V∂σ = zero(eltype(y))
        @simd for kl in 1:IJ
            V∂σ += bs∂σ[kl, ij] * y[kl]
        end
        V∂2μ = zero(eltype(y))
        if (μi != 1) & (μi != μn)
            @simd for kl in 1:IJ
                V∂2μ += bs∂2μ[kl, ij] * y[kl]
            end
        end
        V∂2σ = zero(eltype(y))
        if (σi != 1) & (σi != σn)
            @simd for kl in 1:IJ
                V∂2σ += bs∂2σ[kl, ij] * y[kl]
            end
        end
        if V < 0.0
            ydot[ij] = convert(eltype(y), 1e10)
        else
            ydot[ij] = (byp.ρ * byp.θ * V^(1-1/byp.θ)
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
    @show out.zero
    V = bs' * out.zero
    return μs, σs, V
end



function solve(::Type{Val{:spectral}}, byp::BansalYaronProblem; iterations = 50, kwargs...)
    μmin = byp.μ - 6 * byp.νμ / sqrt(2*byp.κμ)
    μmax = byp.μ + 6 * byp.νμ / sqrt(2*byp.κμ)
    σmin = 1.0 - 3 * byp.νσ / sqrt(2*byp.κσ)
    σmax = 1.0 + 3 * byp.νσ / sqrt(2*byp.κσ)
    if σmin < 0
        σmin = 1e-3
        σmax = 2.0
    end
    basisμ = Basis(Cheb, 20, μmin, μmax)
    basisσ = Basis(Cheb, 20, σmin, σmax)
    basis = Basis(basisμ, basisσ)
    S, (μs, σs) = nodes(basis)
    bs = BasisStructure(basis, Expanded(), S, [0 0 ; 1 0 ; 0 1 ; 2 0 ; 0 2]) 
    bs, bs∂μ, bs∂σ, bs∂2μ, bs∂2σ = bs.vals[1]', bs.vals[2]', bs.vals[3]', bs.vals[4]', bs.vals[5]'
    oldV = bs' \ fill(byp.V, size(bs, 1))
    solve(Val{:spectral}, byp, μs, σs, bs, bs∂μ, bs∂σ, bs∂2μ, bs∂2σ, oldV; iterations = iterations, kwargs...)
end

