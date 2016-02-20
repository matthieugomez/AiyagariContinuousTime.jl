using HJBFiniteDifference, ForwardDiff, Gadfly


##############################################################################
##
## F
##
##############################################################################

function residual(ap::AiyagariProblem, ρπ, σπ, V, g, K, r, w, π, Vdot, gdot, Kdot, rdot, wdot, πdot, Verrors, ηπ)
    newg = vcat(one(eltype(g)) - sum(g), g)
    as = AiyagariSolution(V, newg, K, r, w)
    aa = AiyagariArrays(ap, as)
    HJBFiniteDifference.update_Au!(ap, aa, as)
    HJBFiniteDifference.update_B!(ap, aa, as)
    hjbResidual = aa.u + aa.A' * V - ap.ρ * V .+ Vdot .+ Verrors
    newgIntermediate = aa.A * as.g
    gResidual = gdot -  newgIntermediate[2:end]
    KResidual = K - HJBFiniteDifference.sum_capital(ap.a, newg)
    rResidual = exp(π) * ap.α * K^(ap.α-1) - ap.δ - r
    wResidual = exp(π) * (1-ap.α) * K^ap.α  - w
    πResidual = πdot + (1-ρπ) * π - σπ * ηπ
    return hjbResidual, gResidual, KResidual, rResidual, wResidual, πResidual
end


function unpack(ap::AiyagariProblem, x)
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
    Vdot = x[start:(start + lV - 1)]
    start += lV
    gdot = x[start:(start + lg - 1)]
    start += lg
    Kdot = x[start]
    start += 1
    rdot = x[start]
    start += 1
    wdot = x[start]
    start += 1
    πdot = x[start]
    start += 1
    Verrors = x[start:(start + lV - 1)]
    start += lV
    ηπ = x[start]
    @assert start == length(x)
    return V, g, K, r, w, π, Vdot, gdot, Kdot, rdot, wdot, πdot, Verrors, ηπ
end

function residual(ap, ρπ, σπ, x)
    out = residual(ap, ρπ, σπ, unpack(ap, x)...)
    return vcat(out...)
end

function solve(ap::AiyagariProblem, ρπ, σπ)
    @sprintf "Compute Steady State"
    as = solve(ap)
    g = as.g[2:end]
    ap.invΔ = 0.0
    x0 = vcat(as.V, g, as.K, as.r, as.w, ap.π, zeros(as.V), zeros(g), 0.0, 0.0, 0.0, 0.0, zeros(as.V), 0.0)
    # check function gives 0.0 at steady state 
    @assert maxabs(residual(ap, ρπ, σπ, x0)) < 1e-4

    # compute jacobian around 0
    @sprintf "Compute Jacobian"
    out = jacobian(x -> residual(ap, ρπ, σπ, x), x0)
    # unpack derivatives
    mVarsDerivs = out[:, 1:(length(as.V) + length(g) + 4)]
    mVarsDotDerivs = out[:, (length(as.V) + length(g) + 4 + 1):(2 * length(as.V) + 2 * length(g) + 2 * 4)]
    mEErrorsDerivs = out[:, (2 * length(as.V) + 2 * length(g) + 2 * 4 + 1):(2 * length(as.V) + 2 * length(g) + 2 * 4 + length(as.V))]
    mShocksDerivs = out[:, end:end]

    # solve system
    @sprintf "Solve Linearized Model"
    g0 = mVarsDotDerivs
    g1 = - mVarsDerivs
    c = zeros(length(as.V) + length(g) + 4, 1)
    psi = - mShocksDerivs
    pi = - mEErrorsDerivs
    G1, C, impact, eu = phact_solver(g0, g1, c, psi, pi)
end