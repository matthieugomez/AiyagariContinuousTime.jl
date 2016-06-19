function computeA(x, μ, σ)
    n = length(x)
    Δxm = zeros(x)
    for i in 2:n
        Δxm[i] = x[i] - x[i-1]
    end
    Δxp = zeros(x)
    for i in 1:(n-1)
        Δxp[i] = x[i+1] - x[i]
    end
    Δx = zeros(x)
    for i in 1:n
        Δx[i] = 0.5 * (Δxm[i] + Δxp[i])
    end
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

function kolmogorovforward(A)
    n = size(A, 1)
    # numerical fix 
    for i in 1:n
        A[1, i] = - 1.0
    end
    b = vcat(1.0, zeros(n-1))
    g = - A \ b
    return g
end


function kolmogorovforward(A, δ, ψ)
    n = size(A, 1)
    g = (δ * eye(n) - A) \ (δ * ψ)
    return g
end

kolmogorovforward(x, μ, σ) = kolmogorovforward(computeA(x, μ, σ))
kolmogorovforward(x, μ, σ, args...) = kolmogorovforward(computeA(x, μ, σ)..., args...)