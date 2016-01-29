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
        (ones(μn*σn), ones(μn*σn-1), ones(μn*σn-1), ones((σn-1)*μn), ones((σn-1)*μn)),
        (0, 1, -1, μn, - μn)
    )
    Crows = rowvals(C) 
    Cvals = nonzeros(C) 
    fill!(Cvals, zero(Float64))

    # fill C
    ij = zero(Int)
   @inbounds for σi in 1:σn, μi in 1:μn
        ij += 1
        krange = nzrange(C, ij)
        rows = Crows[krange]
        ∂μ = byp.κμ * (byp.μ - byp.μs[μi]) * byp.invdμ
        ∂2μ = 0.5 * byp.νμ^2 * byp.σs[σi] * byp.invdμ^2
        current =  - min(∂μ, 0.0) + ∂2μ
        if μi > 1
            index = searchsortedfirst(rows, ij - 1)
            Cvals[krange[index]] += current
            index = searchsortedfirst(rows, ij)
            Cvals[krange[index]] -= current
        end

        current =  max(∂μ, 0.0) + ∂2μ
        if μi < μn
            index = searchsortedfirst(rows, ij + 1)
            Cvals[krange[index]] += current
            index = searchsortedfirst(rows, ij)
            Cvals[krange[index]] -= current
        end
        ∂σ = byp.κσ * (1.0 - byp.σs[σi]) * byp.invdσ
        ∂2σ = 0.5 * byp.νσ^2 * byp.σs[σi] * byp.invdσ^2
        current = - min(∂σ, 0.0) + ∂2σ
        if σi > 1
            index = searchsortedfirst(rows, ij - μn)
            Cvals[krange[index]] += current
            index = searchsortedfirst(rows, ij)
            Cvals[krange[index]] -= current
        end

        current = max(∂σ, 0.0) + ∂2σ
        if σi < σn
            index = searchsortedfirst(rows, ij + μn)
            Cvals[krange[index]] += current
            index = searchsortedfirst(rows, ij)
            Cvals[krange[index]] -= current
        end

        current = - byp.ρ * byp.θ + (1-byp.γ) * byp.μs[μi] - 0.5 * (1-byp.γ) * byp.γ * byp.νD^2 * byp.σs[σi]
        index = searchsortedfirst(rows, ij)
        Cvals[krange[index]] += current
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
        (t, y) -> (out = deepcopy(y) ;  f!(byp, C, y, out) ; return out), 
        byp.V, 
        [0.0; 20 * logspace(-1., 7., 9)]; 
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