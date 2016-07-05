function F!{T}(byp::BansalYaronProblem, μs::Vector{Float64}, σs::Vector{Float64}, bs::Matrix{Float64}, bs∂μ::Matrix{Float64}, bs∂σ::Matrix{Float64}, bs∂2μ::Matrix{Float64}, bs∂2σ::Matrix{Float64}, y::Vector{T}, ydot::Vector{T})
    μn = length(μs)
    σn = length(σs)
    IJ = μn * σn
    ij = 0
    @inbounds for σi in 1:σn, μi in 1:μn
        ij += 1
        G = zero(T)
        @simd for kl in 1:IJ
            G += bs[kl, ij] * y[kl]
        end
        G∂μ = zero(T)
        @simd for kl in 1:IJ
            G∂μ += bs∂μ[kl, ij] * y[kl]
        end
        G∂σ = zero(T)
        @simd for kl in 1:IJ
            G∂σ += bs∂σ[kl, ij] * y[kl]
        end
        G∂2μ = zero(T)
        @simd for kl in 1:IJ
            G∂2μ += bs∂2μ[kl, ij] * y[kl]
        end
        G∂2σ = zero(T)
        @simd for kl in 1:IJ
            G∂2σ += bs∂2σ[kl, ij] * y[kl] 
        end
        if (μi == 1 | μi == μn)
            ydot[ij] = G∂μ
        elseif (σi == 1 | σi == σn)
            ydot[ij] = G∂σ
        else
            ydot[ij] = (byp.θ / G - byp.θ * byp.ρ + (1 - byp.γ) * (μs[μi] - 0.5 * byp.γ * byp.νD^2 * σs[σi])
                  + byp.θ  * (byp.κμ * (byp.μ - μs[μi])) * G∂μ / G 
                  + byp.θ  * (byp.κσ * (1.0 - σs[σi])) * G∂σ / G 
                  + 0.5 * byp.θ * byp.νμ^2 * σs[σi] * (G∂2μ / G + (byp.θ - 1) * G∂μ^2 / G^2)
                  + 0.5 * byp.θ * byp.νσ^2 * (G∂2σ / G + (byp.θ - 1) * G∂σ^2 / G^2)
                  )
        end
    end
    return ydot
end

function solve(::Type{Val{:spectral}}, byp::BansalYaronProblem, μs, σs,  bs::Matrix{Float64}, bs∂μ::Matrix{Float64}, bs∂σ::Matrix{Float64}, bs∂2μ::Matrix{Float64}, bs∂2σ::Matrix{Float64}, oldG; maxiter = 100, tol = 1e-10, kwargs...)
    out = nlsolve((y, ydot) -> F!(byp, μs, σs,  bs, bs∂μ, bs∂σ, bs∂2μ, bs∂2σ, y, ydot), oldG; method = :trust_region, show_trace = true, iterations = maxiter, ftol = tol, kwargs...)
    G = bs' * out.zero
    return μs, σs, G
end



function solve(::Type{Val{:spectral}}, byp::BansalYaronProblem; maxiter = 100, tol = 1e-10, kwargs...)
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
    oldG = bs' \ fill(byp.G, size(bs, 1))
    solve(Val{:spectral}, byp, μs, σs, bs, bs∂μ, bs∂σ, bs∂2μ, bs∂2σ, oldG; maxiter = maxiter, tol = tol, kwargs...)
end

