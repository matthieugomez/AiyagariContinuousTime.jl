##############################################################################
##
## Type
##
##############################################################################

type BansalYaronProblemNewton <: BansalYaronProblem

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

    # state grid    
    μs::Vector{Float64}
    σs::Vector{Float64}

    # arrays for storage
    C::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}
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

function BansalYaronProblemNewton(;μ = 0.0015, νD = 0.0078, κμ = 0.0212, κσ = 0.0131, νμ = 0.0003432, νσ = 0.0378,  ρ = 0.002, γ = 7.5, ψ = 1.5, μn = 100, σn = 100)

    θ = (1-γ)/(1-1/ψ)

    # create drift grid +-6 sd of stationary distribution
    μmin = μ - 6 * νμ / sqrt(2*κμ)
    μmax = μ + 6 * νμ / sqrt(2*κμ)
    μs = collect(linspace(μmin, μmax, μn))
    invdμ = (μn - 1)/(μmax - μmin)

    # create volatility +-3 sd. if σmin <0, just do 1e-3, 2.0
    σmin = 1.0 - 3 * νσ / sqrt(2*κσ)
    σmax = 1.0 + 3 * νσ / sqrt(2*κσ)
    if σmin < 0
        σmin = 1e-3
        σmax = 2.0
    end
    σs = collect(linspace(σmin, σmax, σn))
    invdσ = (σn - 1)/(σmax - σmin)

    # create sparse matrix of the right form finite difference matrix B 
    C = spdiagm(
        (ones(μn*σn), ones(μn*σn-1), ones(μn*σn-1), ones((σn-1)*μn), ones((σn-1)*μn)),
        (0, 1, -1, μn, - μn)
    )
    Cvals = nonzeros(C) 
    fill!(Cvals, zero(Float64))

    # set C such that C'v^n+1 = u^n
    Crows = rowvals(C) 
    ij = zero(Int)
    @inbounds for σi in 1:σn, μi in 1:μn
        ij += 1
        krange = nzrange(C, ij)
        rows = Crows[krange]
        current =  min(κμ * (μ - μs[μi]) * invdμ, 0.0) - 0.5 * νμ^2 * σs[σi] * invdμ^2
        if μi > 1
            index = searchsortedfirst(rows, ij - 1)
            Cvals[krange[index]] += current
            index = searchsortedfirst(rows, ij)
            Cvals[krange[index]] -= current
        end

        current =  - max(κμ * (μ - μs[μi]) * invdμ, 0.0) - 0.5 * νμ^2 * σs[σi] * invdμ^2
        if μi < μn
            index = searchsortedfirst(rows, ij + 1)
            Cvals[krange[index]] += current
            index = searchsortedfirst(rows, ij)
            Cvals[krange[index]] -= current
        end

        current = min(κσ * (1.0 - σs[σi]) * invdσ, 0.0) - 0.5 * νσ^2 * σs[σi] * invdσ^2
        if σi > 1
            index = searchsortedfirst(rows, ij - μn)
            Cvals[krange[index]] += current
            index = searchsortedfirst(rows, ij)
            Cvals[krange[index]] -= current
        end

        current = - max(κσ * (1.0 - σs[σi]) * invdσ, 0.0) - 0.5 * νσ^2 * σs[σi] * invdσ^2
        if σi < σn
            index = searchsortedfirst(rows, ij + μn)
            Cvals[krange[index]] += current
            index = searchsortedfirst(rows, ij)
            Cvals[krange[index]] -= current
        end
        
        current = ρ * θ - (1-γ) * μs[μi] + 0.5 * (1-γ) * γ * νD^2 * σs[σi]
        index = searchsortedfirst(rows, ij)
        Cvals[krange[index]] += current
    end
    B = deepcopy(C)
    # initialize value at stationary value
    V = Array(Float64, μn * σn)
    fill!(V, (-1/(θ*ρ) * (μ * (1-γ) - 0.5 * (1-γ) * γ * νD^2 * 1.0) + 1.0)^(-1/(1-1/θ)))
    newV = deepcopy(V)
    u = fill(zero(Float64), μn*σn)
    BansalYaronProblemNewton(μ, νD, κμ , κσ , νμ , νσ, ρ, γ, ψ, θ, μs, σs, C, B, V, newV, u)
end

##############################################################################
##
## Solve
##
##############################################################################

function update_value!(byp::BansalYaronProblemNewton)
    μn = length(byp.μs)
    σn = length(byp.σs)

    # update B
    copy!(nonzeros(byp.B), nonzeros(byp.C))
    Brows = rowvals(byp.B)
    Bvals = nonzeros(byp.B)
    ij = zero(Int)
    @inbounds for σi in 1:σn, μi in 1:μn
        ij += 1
        krange = nzrange(byp.B, ij)
        rows = Brows[krange]
        # so that sum of C by column equal zero
        current = - byp.ρ * (byp.θ-1) * byp.V[ij]^(-1/byp.θ)
        index = searchsortedfirst(rows, ij)
        Bvals[krange[index]] += current
    end


    # update U
    ij = zero(Int)
    @inbounds for σi in 1:σn, μi in 1:μn
        ij += 1
        byp.u[ij] = byp.ρ * byp.V[ij]^(1-1/byp.θ)
    end

    # solve for new V
    byp.newV = byp.B' \ byp.u
end



