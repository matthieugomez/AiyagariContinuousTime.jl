function computeA(x::AbstractVector, μ::AbstractVector, σ::AbstractVector)
    n = length(x)
    # construct small Δx for finite difference scheme
    Δxm = zeros(n)
    for i in 2:n
        Δxm[i] = x[i] - x[i-1]
    end
    Δxp = zeros(n)
    for i in 1:(n-1)
        Δxp[i] = x[i+1] - x[i]
    end
    Δx = zeros(n)
    for i in 1:n
        Δx[i] = 0.5 * (Δxm[i] + Δxp[i])
    end
    # construct matrix A. The key is that sum of each column = 0.0 and off diagonals are positive (singular M-matrix)
    A = zeros(n, n)
    for i in 1:length(x)
        # drift
        if (μ[i] < 0) & (i > 1)
            A[i - 1, i] +=  - μ[i] / Δxm[i]
            A[i, i] +=  μ[i] / Δxm[i]
        elseif (μ[i] >= 0) & (i < n)
            A[i, i] +=  - μ[i] / Δxp[i]
            A[i + 1, i] +=  μ[i] / Δxp[i]
        end
        # volatility
        if (i > 1)
            A[i-1, i] += 0.5 * σ[i]^2 / (Δx[i] * Δxm[i])
            A[i, i] += - 0.5 * σ[i]^2 / (Δx[i] * Δxm[i])
        end
        if (i < n)
            A[i+1, i] += 0.5 * σ[i]^2 / (Δx[i] * Δxp[i])
            A[i, i] += - 0.5 * σ[i]^2 / (Δx[i] * Δxp[i])
        end
    end
    return A
end

function kolmogorovforward(A::Matrix)
    # transformation is like adding sum = 1.0
    n = size(A, 2)
    for j in 1:n
        A[1, j] = 1e-10
    end
    b = vcat(1e-10, zeros(n-1))
    density = A \ b
    # some values are -eps() < density < 0, — just set them to zero
    density = abs(density) ./ sumabs(density)
end


function kolmogorovforward(A::Matrix, δ::Real, ψ::AbstractVector)
    n = size(A, 1)
    (δ * eye(n) - A) \ (δ * ψ)
end

kolmogorovforward(x::AbstractVector, μ::AbstractVector, σ::AbstractVector) = kolmogorovforward(computeA(x, μ, σ))
kolmogorovforward(x::AbstractVector, μ::AbstractVector, σ::AbstractVector, args...) = kolmogorovforward(computeA(x, μ, σ), args...)
