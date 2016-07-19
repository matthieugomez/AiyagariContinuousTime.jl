function make_grid(x::AbstractVector)
    n = length(x)
    # construct small Δx for finite difference scheme
    Δxm = zeros(x)
    for i in 2:n
        Δxm[i] = x[i] - x[i-1]
    end
    Δxp = zeros(x)
    for i in 1:(n-1)
        Δxp[i] = x[i+1] - x[i]
    end
    Δx = zeros(x)
    # difference with Moll notes so that gives usual one when grid is uniform
    Δx[1] = Δxp[1]
    Δx[end] = Δxm[end]
    for i in 2:(n-1)
        Δx[i] = 0.5 * (Δxm[i] + Δxp[i])
    end
    return Δx, Δxm, Δxp
end

function computeA(x::AbstractVector, μ::AbstractVector, σ::AbstractVector)
    n = length(x)
    T = eltype(μ)
    Δx, Δxm, Δxp = make_grid(x)
    # construct matrix A. The key is that sum of each column = 0.0 and off diagonals are positive (singular M-matrix)
    A = zeros(T, n, n)
    for i in 1:n
        # drift
        # if drift positive at the boundary, we use downard derivative (so that mean is conserved in stationary distribution => important to satisfy walras law)
        # none is due to non unifrm grid
        if ((μ[i] < 0) & (i > 1)) | (i == n)
            current = μ[i] / Δxm[i]
            A[i-1, i] -=  current
            A[i, i] +=  current
        elseif (μ[i] >= 0) & (i < n)
            current = μ[i] / Δxp[i]
            A[i, i] -= current
            A[i+1, i] += current
        end
        # volatility
        if (i > 1) && (i < n)
            current = 0.5 * σ[i]^2 / (Δx[i] * Δxm[i])
            A[i-1, i] += current
            A[i, i] -= current
   
            current = 0.5 * σ[i]^2 / (Δx[i] * Δxp[i])
            A[i+1, i] += current
            A[i, i] -= current
        end
    end
    #check column sum to zero
    for i in 1:n
        @assert abs(sum(slice(A, :, i))) <= 1e-10
        #important for walras law
        #two issues : μ positive at the top
        #@show (sum(slice(A, :, i) .* x) - μ[i])/(max(1e-5, μ[i]))
    end
    @show 
    return A
end

function clean!(g)
    # in case some values are -eps() < g < 0, set them to zero
    @assert all(g .> -1e-5)
    map!(x -> max(x, zero(eltype(g))), g)
    scale!(g,  1 / sum(g))
end


function kolmogorovforward(A::AbstractMatrix)
    n = size(A, 2)
    for j in 1:n
        A[1, j] = 1.0
    end
    b = vcat(1.0, zeros(n-1))
    g = A \ B
    clean!(g)
end


function kolmogorovforward(A::AbstractMatrix, δ::Union{AbstractVector, Real}, ϕ::AbstractVector)
    n = size(A, 1)
    g = (δ .* eye(n) .- A) \ (δ .* ϕ)
    A, clean!(g)
end

kolmogorovforward(x::AbstractVector, μ::AbstractVector, σ::AbstractVector) = kolmogorovforward(computeA(x, μ, σ))
kolmogorovforward(x::AbstractVector, μ::AbstractVector, σ::AbstractVector, args...) = kolmogorovforward(computeA(x, μ, σ), args...)