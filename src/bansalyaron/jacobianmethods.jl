##############################################################################
##
## Structure
## F_i(y) can be written C x y + f(y_{i}) = 0
## The Jacobian is simply C' + diag(f')
##############################################################################

function structure(byp::BansalYaronProblem)
    μn = length(byp.μs)
    σn = length(byp.σs)
    # set up C
    C = spdiagm(
        (ones(μn*σn), ones(μn*σn-1), ones(μn*σn-1), ones(μn*σn-2), ones(μn*σn-2), ones((σn-1)*μn), ones((σn-1)*μn), ones((σn-2)*μn), ones((σn-2)*μn)),
        (0, 1, -1, 2, -2, μn, - μn, 2 * μn, - 2 * μn)
    )
    fill!(nonzeros(C), zero(Float64))

    # fill C
    ij = zero(Int)
    @inbounds for σi in 1:σn, μi in 1:μn
        ij += 1
        ∂μ = byp.κμ * (byp.μ - byp.μs[μi]) * byp.invdμ
        ∂2μ = 0.5 * byp.νμ^2 * byp.σs[σi] * byp.invdμ^2
        if μi == 1
            C[ij, ij] += - 0.5 * 3 * ∂μ
            C[ij + 1, ij] +=  0.5 * 4 * ∂μ
            C[ij + 2, ij] += - 0.5 * ∂μ
        elseif μi == μn
            C[ij, ij] += 0.5 * 3 * ∂μ
            C[ij - 1, ij] += - 0.5 * 4 * ∂μ
            C[ij - 2, ij] += 0.5 * ∂μ
        else
            C[ij - 1, ij] += - ∂μ * 0.5 + ∂2μ
            C[ij + 1, ij] += ∂μ * 0.5 + ∂2μ
            C[ij, ij] += - 2 * ∂2μ
        end
        ∂σ = byp.κσ * (1.0 - byp.σs[σi]) * byp.invdσ
        ∂2σ = 0.5 * byp.νσ^2 * byp.σs[σi] * byp.invdσ^2
        if σi == 1
            C[ij, ij] +=  - 0.5 * 3  * ∂σ
            C[ij + μn, ij] += 0.5 * 4 * ∂σ
            C[ij + 2 * μn, ij] += - 0.5 * ∂σ
        elseif σi == σn
            C[ij, ij] += 0.5 * 3 * ∂σ
            C[ij - μn, ij] += - 0.5 * 4 * ∂σ
            C[ij - 2 * μn, ij] += 0.5 * ∂σ
        else
            C[ij - μn, ij] += - ∂σ * 0.5 + ∂2σ
            C[ij + μn, ij] += ∂σ * 0.5 + ∂2σ
            C[ij, ij] += - 2 * ∂2σ
        end
        C[ij, ij] += - byp.ρ * byp.θ + (1-byp.γ) * byp.μs[μi] - 0.5 * (1-byp.γ) * byp.γ * byp.νD^2 * byp.σs[σi]
    end
    return C
end

# Express F(y) from C
function f!(byp::BansalYaronProblem, C, x::Vector{Float64}, out::Vector{Float64})
    ij = zero(Int)
    μn = length(byp.μs)
    σn = length(byp.σs)
    Cvals = nonzeros(C)
    Crows = rowvals(C)
    @inbounds for σi in 1:σn, μi in 1:μn
        ij += 1
        krange = nzrange(C, ij)
        rows = Crows[krange]
        current = zero(Float64)
        for k in krange
            current += Cvals[k] * x[Crows[k]]
        end
        current += byp.ρ * byp.θ * max(x[ij], 0.0)^(1-1/byp.θ)
        out[ij] = current
    end
end

function f(byp::BansalYaronProblem, C, x::Vector{Float64})
    out = deepcopy(x)
    f!(byp, C, x, out)
    return out
end

# Express Jac(F)(y) from Ct
function g!(byp::BansalYaronProblem, Ct, x::Vector{Float64}, out)
    ij = zero(Int)
    μn = length(byp.μs)
    σn = length(byp.σs)
    Ctvals = nonzeros(Ct)
    outvals = nonzeros(out)
    outrows = rowvals(out)
    ij = zero(Int)
    @inbounds for σi in 1:σn, μi in 1:μn
        ij += 1
        for k in nzrange(Ct, ij)
            outvals[k] = Ctvals[k]
            if outrows[k] == ij
                outvals[k] +=  byp.ρ * (byp.θ -1) * max(x[ij], 0.0)^(-1/byp.θ)
            end
        end
    end
    return out
end

function g(byp::BansalYaronProblem, Ct, x::Vector{Float64})
    out = deepcopy(Ct)
    g!(byp, Ct, x, out)
    return out
end

##############################################################################
##
## Solution methods
##
##############################################################################

function solve(::Type{Val{:ode23s}}, byp::BansalYaronProblem)
    C = structure(byp)
    Ct = C'
    tout, yout = ODE.ode23s(
        (t, y) -> f(byp, C, y),
        byp.V, 
        [0.0; 10.0]; 
        jacobian = (t, y) -> g(byp, Ct, y))
    return yout[end]
end

function solve(::Type{Val{:newton}}, byp::BansalYaronProblem)
    C = structure(byp)
    Ct = C'
    oldV = deepcopy(byp.V)
    u = similar(oldV)
    jac = deepcopy(Ct)
    for iter in 1:1_000
        g!(byp, Ct, oldV, jac)
        for i in eachindex(oldV)
            u[i] = byp.ρ * oldV[i]^(1-1/byp.θ)
        end
        newV = jac \ u
        scale!(newV, -1)
        distance = chebyshev(vec(newV), vec(oldV))
        @show iter, distance
        if distance < 1e-10
            return newV
        else
            oldV, newV = newV, oldV
        end
    end
end

function solve(::Type{Val{:trust_region}}, byp::BansalYaronProblem)
    C = structure(byp)
    Ct = C'
    df = DifferentiableGivenSparseMultivariateFunction(
        (x, out) -> f!(byp, C, x, out), 
        (x, out)-> g!(byp, Ct, x, out),
        Ct
        )
    out = nlsolve(df, byp.V, method = :trust_region, show_trace = true, xtol = 1e-6, ftol = 1e-15)
    return out.zero
end