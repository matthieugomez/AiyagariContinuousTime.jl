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

function computeA(x::AbstractVector, μ::AbstractVector, σ::AbstractVector; check = true)
    n = length(x)
    T = eltype(μ)
    Δx, Δxm, Δxp = make_grid(x)
    # construct matrix A. The key is that sum of each column = 0.0 and off diagonals are positive (singular M-matrix)
    A = zeros(T, n, n)
    for i in 1:n
        # drift
        # if drift positive at the boundary, we use downard derivative (so that mean is conserved in stationary distribution => important to satisfy walras law)
        # none is due to non unifrm grid
        if ((μ[i] < 0) && (i > 1)) || (i == n)
            current = μ[i] / Δxm[i]
            A[i-1, i] -=  current
            A[i, i] +=  current
        elseif ((μ[i] >= 0) && (i < n)) || (i == 1)
            current = μ[i] / Δxp[i]
            A[i, i] -= current
            A[i+1, i] += current
        end
        # volatility
        if (i > 1) && (i <n)
            current = 0.5 * σ[i]^2 / (Δx[i] * Δxm[i])
            A[i-1, i] += current
            A[i, i] -= current
    
            current = 0.5 * σ[i]^2 / (Δx[i] * Δxp[i])
            A[i+1, i] += current
            A[i, i] -= current
        end
    end
    #check column sum to zero
    if check
        v1 = vec(sum(A, 1))
        if any(abs(v1) .> 1e-10)
            info("columns don't sum up to one")
        end
        v2 = vec(sum(A .* x, 1)) ./ μ .- 1.0
        if any(abs(v2) .> 1e-10)
            info("aggregate wealth not conserved")
        end
    end
    return A
end

function clean!(g)
    # in case some values are -eps() < g < 0, set them to zero
    if any(g .< 0.0)
        @assert all(g .> -1e-5)
        map!(x -> max(x, zero(eltype(g))), g)
        scale!(g,  1 / sum(g))
    end
    return g
end



# Stationary
function stationary(A::AbstractMatrix; check = true)
    n = size(A, 2)
    for j in 1:n
        A[1, j] = 1.0
    end
    b = vcat(1.0, zeros(n-1))
    g = A \ B
    if check
        clean!(g)
    end
    return g
end
function stationary(A::AbstractMatrix, δ::Union{AbstractVector, Real}, ϕ::AbstractVector; check = true)
    n = size(A, 2)
    g = (δ .* eye(n) .- A) \ (δ .* ϕ)
    if check
        clean!(g)
    end
    return g
end
stationary(x::AbstractVector, μ::AbstractVector, σ::AbstractVector; check = true) = stationary(computeA(x, μ, σ; check = check); check = check)
stationary(x::AbstractVector, μ::AbstractVector, σ::AbstractVector, args...; check = true) = stationary(computeA(x, μ, σ, check = check), args...; check = check)


# forward
function forward(g, A::AbstractMatrix, δ::Union{AbstractVector, Real}, ϕ::AbstractVector)
    A * g  - δ .* g + δ .* ϕ
end
function forward(g, A::AbstractMatrix)
    A * g
end

function backward(g, A::AbstractMatrix, δ::Union{AbstractVector, Real}, ϕ::AbstractVector)
    A * g  - δ .* g + δ .* ϕ
end
function backward(g, A::AbstractMatrix)
    A * g
end


forward(g::AbstractVector, x::AbstractVector, μ::AbstractVector, σ::AbstractVector) = forward(g, computeA(x, μ, σ))
forward(g::AbstractVector, x::AbstractVector, μ::AbstractVector, σ::AbstractVector, check = true) = forward(g, computeA(x, μ, σ), check = check)
