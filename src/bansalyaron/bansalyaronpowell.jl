##############################################################################
##
## Type
##
##############################################################################

type BansalYaronProblemPowell <: BansalYaronProblem

    # consumption process parameters
    μ::Float64 
    νD::Float64
    κμ::Float64 
    κσ::Float64 
    νμ::Float64 
    νσ::Float64 

    # utility parameters
    ρ::Float64  
    γ::Float64 
    ψ::Float64
    θ::Float64

    # algorithm parameters
    invdμ::Float64
    invdσ::Float64 

    # grid    
    μs::Vector{Float64}
    σs::Vector{Float64}

    # storage array
    C::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}
    Ct::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}
    B::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}

    V::Vector{Float64}
    newV::Vector{Float64}    
    u::Vector{Float64}
end


##############################################################################
##
## Constructor
##
##############################################################################


function BansalYaronProblemPowell(;μ = 0.0015, νD = 0.0078, κμ = 0.0212, κσ = 0.0131, νμ = 0.0003432, νσ = 0.0378,  ρ = 0.002, γ = 7.5, ψ = 1.5, μn = 100, σn = 100)

    θ = (1-γ)/(1-1/ψ)

    # create drift grid +-3 sd of stationary distribution
    μmin = μ - 6 * νμ / sqrt(2*κμ)
    μmax = μ + 6 * νμ / sqrt(2*κμ)
    μs = collect(linspace(μmin, μmax, μn))
    invdμ = (μn - 1)/(μmax - μmin)

    # create volatility +-3 sd otherwise negative volatility
    σmin = 1.0 - 3 * νσ / sqrt(2*κσ)
    σmax = 1.0 + 3 * νσ / sqrt(2*κσ)
    if σmin < 0
        σmin = 1e-3
        σmax = 2.0
    end
    σs = collect(linspace(σmin, σmax, σn))
    invdσ = (σn - 1)/(σmax - σmin)

    # create matrix of the right form filled with zeros
    C = spdiagm(
        (ones(μn*σn), ones(μn*σn-1), ones(μn*σn-1), ones((σn-1)*μn), ones((σn-1)*μn)),
        (0, 1, -1, μn, - μn)
    )

    # set up C
    Crows = rowvals(C) 
    Cvals = nonzeros(C) 
    fill!(Cvals, zero(Float64))
    ij = zero(Int)
    @inbounds for σi in 1:σn, μi in 1:μn
        ij += 1
        krange = nzrange(C, ij)
        rows = Crows[krange]
        ∂μ = κμ * (μ - μs[μi]) * invdμ
        ∂2μ = 0.5 * νμ^2 * σs[σi] * invdμ^2
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
        ∂σ = κσ * (1.0 - σs[σi]) * invdσ
        ∂2σ = 0.5 * νσ^2 * σs[σi] * invdσ^2
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
        
        current = - ρ * θ + (1-γ) * μs[μi] - 0.5 * (1-γ) * γ * νD^2 * σs[σi]
        index = searchsortedfirst(rows, ij)
        Cvals[krange[index]] += current
    end
    Ct = C'
    B = deepcopy(Ct)
    # initialize value at stationary value
    V = Array(Float64, μn * σn)
    fill!(V, (-1/(θ*ρ) * (μ * (1-γ) - 0.5 * (1-γ) * γ * νD^2 * 1.0) + 1.0)^(-1/(1-1/θ)))
    newV = deepcopy(V)
    u = fill(zero(Float64), μn*σn)
    BansalYaronProblemPowell(μ, νD, κμ , κσ , νμ , νσ, ρ, γ, ψ, θ, invdμ, invdσ, μs, σs, C, Ct, B, V, newV, u)
end

##############################################################################
##
## Solve
##
##############################################################################


function f!(byp::BansalYaronProblemPowell, x::Vector{Float64}, out::Vector{Float64})
    ij = zero(Int)
    μn = length(byp.μs)
    σn = length(byp.σs)
    Cvals = nonzeros(byp.C)
    Crows = rowvals(byp.C)
    @inbounds for σi in 1:σn, μi in 1:μn
        ij += 1
        krange = nzrange(byp.C, ij)
        rows = Crows[krange]
        current = zero(Float64)
        for k in krange
            current += Cvals[k] * x[Crows[k]]
        end
        current += byp.ρ * byp.θ * max(x[ij], 0.0)^(1-1/byp.θ)
        out[ij] = current
    end
end


function g!(byp::BansalYaronProblemPowell, x::Vector{Float64}, out::Base.SparseMatrix.SparseMatrixCSC{Float64, Int})
    ij = zero(Int)
    μn = length(byp.μs)
    σn = length(byp.σs)
    Ctvals = nonzeros(byp.Ct)
    Bvals = nonzeros(byp.B)
    Brows = rowvals(byp.B)
    ij = zero(Int)
    @inbounds for σi in 1:σn, μi in 1:μn
        ij += 1
        for k in nzrange(byp.B, ij)
            Bvals[k] = Ctvals[k]
            if Brows[k] == ij
                Bvals[k] +=  byp.ρ * byp.θ * (1-1/byp.θ) * max(x[ij], 0.0)^(-1/byp.θ)
            end
        end
    end
    out.colptr = byp.B.colptr
    out.rowval = byp.B.rowval
    out.nzval = byp.B.nzval
end


function update_value!(byp::BansalYaronProblemPowell; method::Symbol = :trust_region, iterations = 10_000)
    μn = length(byp.μs)
    σn = length(byp.σs)
    df = DifferentiableSparseMultivariateFunction(
        (x, out) -> f!(byp, x, out), 
        (x, out)-> g!(byp, x, out)
        )
    out = nlsolve(df, byp.V, method = method, show_trace = true, iterations = iterations, xtol = 1e-6, ftol = 1e-15)
    @assert all(x -> x >= zero(Float64), out.zero)
    byp.newV = out.zero
end
