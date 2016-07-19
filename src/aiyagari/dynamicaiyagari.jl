##############################################################################
##
## DynamicSolution
##
##############################################################################

type DynamicAiyagariSolution{T}
    as::AiyagariSolution{T}
    G1::Matrix{T}
    impact::Matrix{T}
    existence::Bool
    uniqueness::Bool
end


##############################################################################
##
## F
##
##############################################################################

function F(ap::AiyagariProblem, ρπ, σπ, V, g, K, r, w, π, Vpost, gpost, Kpost, rpost, wpost, πpost, Verrors, ηπ)
    newg = vcat(one(eltype(g)) - sum(g), g)
    as = AiyagariSolution(V, newg, K, r, w)
    aa = AiyagariArrays(ap, as)
    update_Au!(ap, aa, as)
    update_B!(ap, aa, as)
    hjbResidual = aa.u + full(aa.A)' * V - ap.ρ * V .+ Vpost .- V .+ Verrors
    newgIntermediate = full(aa.A) * as.g
    gResidual = gpost .-  g .- newgIntermediate[2:end]
    KResidual = K - HJBFiniteDifference.sum_capital(ap.a, newg)
    rResidual = exp(π) * ap.α * K^(ap.α-1) - ap.δ - r
    wResidual = exp(π) * (1-ap.α) * K^ap.α  - w
    πResidual = πpost .- π .- (1-ρπ) * π - σπ * ηπ
    return hjbResidual, gResidual, KResidual, rResidual, wResidual, πResidual
end


function F(ap, ρπ, σπ, x)
    lV = length(ap.a) * length(ap.z)
    lg = length(ap.a) * length(ap.z) - 1
    start = 1
    V = x[start:(start + lV - 1)]
    start += lV
    g = x[start:(start + lg - 1)]
    start += lg
    K = x[start]
    start += 1
    r = x[start]
    start += 1
    w = x[start]
    start += 1
    π = x[start]
    start += 1
    Vpost = x[start:(start + lV - 1)]
    start += lV
    gpost = x[start:(start + lg - 1)]
    start += lg
    Kpost = x[start]
    start += 1
    rpost = x[start]
    start += 1
    wpost = x[start]
    start += 1
    πpost = x[start]
    start += 1
    Verrors = x[start:(start + lV - 1)]
    start += lV
    ηπ = x[start]
    @assert start == length(x)
    return vcat(F(ap, ρπ, σπ, V, g, K, r, w, π, Vpost, gpost, Kpost, rpost, wpost, πpost, Verrors, ηπ)...)
end

function solve(ap::AiyagariProblem, ρπ, σπ)
    println("Compute Steady State")
    as = solve(ap, verbose = false)
    g = as.g[2:end]
    x0 = vcat(as.V, g, as.K, as.r, as.w, ap.π, deepcopy(as.V), deepcopy(g), as.K, as.r, as.w, ap.π, zeros(as.V), 0.0)
    ap.invΔ = 0.0

    # check function gives 0.0 at steady state 
    @assert maxabs(F(ap, ρπ, σπ, x0)) < 1e-4

    # compute jacobian around 0
    println("Compute Jacobian")
    out = ForwardDiff.jacobian(x -> F(ap, ρπ, σπ, x), x0)
    # unpack derivatives
    Γ1 = - out[:, 1:(length(as.V) + length(g) + 4)]
    Γ0 = out[:, (length(as.V) + length(g) + 4 + 1):(2 * length(as.V) + 2 * length(g) + 2 * 4)]
    Π = - out[:, (2 * length(as.V) + 2 * length(g) + 2 * 4 + 1):(2 * length(as.V) + 2 * length(g) + 2 * 4 + length(as.V))]
    Ψ = - out[:, end:end]

    # solve system
    println("Solve Linearized Model")
    c = zeros(length(as.V) + length(g) + 4, 1)
    G1, C, impact, fmat, fwt, ywt, gev, eu, loose = gensys(Γ0, Γ1, c, Ψ, Π)
    existence = eu[1] == 1
    uniquenes = eu[2] == 1
    return DynamicAiyagariSolution(as, G1, impact, existence, uniqueness)
end


function simulate(ap::AiyagariProblem, as::DynamicAiyagariSolution)
    T = 200
    N = 2000
    T1 = div(N, T)
    dt = T / N
    lV = length(ap.a) * length(ap.z)
    lg = length(ap.a) * length(ap.z) - 1
    g0 = as.as.g[2:end]
    x = zeros(Float64, lV + lg + 4)
    V = Array(Float64, lV, N)
    g = Array(Float64, lg, N)
    K = Array(Float64, N)
    r = Array(Float64, N)
    w = Array(Float64, N)
    A = as.G1 * dt + speye(length(x))
    B = as.impact * sqrt(dt)
    for n in 1:N
        # unpack solution
        start = 1
        V[:, n] = x[start:(start + lV - 1)]
        start += lV
        g[:, n] = x[start:(start + lg - 1)]
        start += lg
        K[n] = x[start]
        start += 1
        r[n] = x[start]
        start += 1
        w[n] = x[start]
        if n <= T1
            newx = A * x + B
        else
            newx = A * x
        end

        # update
        newx, x = x, newx
    end
    return linspace(1, T, N), V, g, K/as.as.K, r/as.as.r, w/as.as.w
end