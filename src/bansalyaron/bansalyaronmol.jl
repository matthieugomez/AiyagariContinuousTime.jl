##############################################################################
##
## Type
##
##############################################################################

type BansalYaronProblemMOL 
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

    V::Vector{Float64}
    newV::Vector{Float64}
    jac::Matrix{Float64}
end


##############################################################################
##
## Constructor
##
##############################################################################


function BansalYaronProblemMOL(;μ = 0.0015, νD = 0.0078, κμ = 0.0212, κσ = 0.0131, νμ = 0.0003432, νσ = 0.0378,  ρ = 0.002, γ = 7.5, ψ = 1.5, μn = 100, σn = 100)

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

    # initialize value at stationary value
    V = Array(Float64, μn * σn)
    fill!(V, (-1/(θ*ρ) * (μ * (1-γ) - 0.5 * (1-γ) * γ * νD^2 * 1.0) + 1.0)^(-1/(1-1/θ)))

    BansalYaronProblemMOL(μ, νD, κμ , κσ , νμ , νσ, ρ, γ, ψ, θ, invdμ, invdσ, μs, σs, V)
end

##############################################################################
##
## Solve
##
##############################################################################
function shift(x::Vector, i, s)
    newi = i + s
    if (newi <= length(x)) & (newi >= 1)
        return x[newi]
    else
        return x[i]
    end
end

function F!(byp::BansalYaronProblemMOL, y::Vector, ydot::Vector)
    ij = zero(Int)
    μn = length(byp.μs)
    σn = length(byp.σs)
    @inbounds for σi in 1:σn, μi in 1:μn
        ij += 1
        ∂μ = byp.κμ * (byp.μ - byp.μs[μi]) * byp.invdμ
        ∂2μ = 0.5 * byp.νμ^2 * byp.σs[σi] * byp.invdμ^2
        ∂σ = byp.κσ * (1.0 - byp.σs[σi]) * byp.invdσ
        ∂2σ = 0.5 * byp.νσ^2 * byp.σs[σi] * byp.invdσ^2
        current = 0.0
        current += ∂2μ * (shift(y, ij, 1) + shift(y, ij, -1) - 2 * shift(y, ij, 0)) 
        if ∂μ > 0 
            current += ∂μ * (shift(y, ij, 1) - shift(y, ij, 0))
        else 
            current += ∂μ * (shift(y, ij, 0) - shift(y, ij, -1))
        end
        current += ∂2σ * (shift(y, ij, μn) + shift(y, ij, -μn) - 2 * shift(y, ij, 0))
        if ∂σ > 0 
          current += ∂σ * (shift(y, ij, μn) - shift(y, ij, 0))
        else 
          current += ∂σ * (shift(y, ij, 0) - shift(y, ij, - μn))
        end 
        current += byp.ρ * byp.θ * max(y[ij], 0.0)^(1-1/byp.θ) + (- byp.ρ * byp.θ + (1-byp.γ) * byp.μs[μi] - 0.5 * (1-byp.γ) * byp.γ * byp.νD^2 * byp.σs[σi])*y[ij] 
        ydot[ij] = current
    end
    @show norm(ydot)
    return ydot
end

function F(byp::BansalYaronProblemMOL, y::Vector)
    ydot = deepcopy(y)
    F!(byp, y, ydot)
end


function solve(byp::BansalYaronProblemMOL; kwargs...)
    oldV = deepcopy(byp.V)
    for iter in 1:10
       tout, yout = ODE.ode23((t, y) -> F(byp, y), oldV, [0.0;100.0]; kwargs...)
       newV = yout[end]
       distance = chebyshev(vec(oldV), vec(newV))
       if distance < 1e-10
            return newV
        else
            @show iter distance
            oldV, newV = newV, oldV
        end
    end
    return oldV
end


#function solve(byp::BansalYaronProblemMOL; kwargs....)
#    res = Sundials.cvode((t, y, ydot) -> F!(byp, y, ydot), byp.V, [0.0;100.0])
#    return vec(yout[end, :])
#end



# using Gadfly, HJBFiniteDifference ; bypp = BansalYaronProblem(ρ = -log(0.9989), γ = 7.5, ψ = 1.5, νD = 0.0072, νμ = 0.038 * 0.0072, νσ = 0.0000028 / 0.0072^2, κμ = -log(0.975), κσ = -log(0.999)) ; solution = solve(bypp) ; plot(bypp, solution, :s2) ; byp = BansalYaronProblemMOL(ρ = -log(0.9989), γ = 7.5, ψ = 1.5, νD = 0.0072, νμ = 0.038 * 0.0072, νσ = 0.0000028 / 0.0072^2, κμ = -log(0.975), κσ = -log(0.999)) ; solve(byp)
# copy!(byp.V, solution)
# solution = solve(byp)
# maxabs(HJBFiniteDifference.F(byp, solution))
