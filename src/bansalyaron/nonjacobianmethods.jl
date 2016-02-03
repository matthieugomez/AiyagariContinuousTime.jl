##############################################################################
##
## Solve
##
##############################################################################
type Grid{T}
    x::Vector{T}
    I::Int
    J::Int
end

function Base.getindex(x::Grid, i, j)
    return x.x[i + x.I * (j - 1)]
end


function F!(byp::BansalYaronProblem, y::Vector, ydot::Vector)
    μn = length(byp.μs)
    σn = length(byp.σs)
    fy = Grid(y, μn, σn)
    ij = 0
    for σi in 1:σn, μi in 1:μn
        ij += 1
        ∂μ = byp.κμ * (byp.μ - byp.μs[μi]) * byp.invdμ
        ∂2μ = 0.5 * byp.νμ^2 * byp.σs[σi] * byp.invdμ^2
        ∂σ = byp.κσ * (1.0 - byp.σs[σi]) * byp.invdσ
        ∂2σ = 0.5 * byp.νσ^2 * byp.σs[σi] * byp.invdσ^2
        v∂σ = 0.0
        v∂2σ = 0.0
        if  (σi > 1) & (σi < σn)
            v∂2σ = (fy[μi, σi+1] + fy[μi, σi-1] - 2.0 * fy[μi, σi]) 
            v∂σ = 0.5 * (fy[μi, σi+1] - fy[μi, σi-1]) 
        elseif σi == 1
            v∂σ = - 0.5 * (3 * fy[μi, σi] - 4 * fy[μi, σi+1] + fy[μi, σi+2])
        elseif σi == σn
            v∂σ = 0.5 * (3 * fy[μi, σi] - 4 * fy[μi, σi-1] + fy[μi, σi-2])
        end
        v∂μ = 0.0
        v∂2μ = 0.0
        if  (μi > 1) & (μi < μn)
            v∂2μ = (fy[μi+1, σi] + fy[μi-1, σi] - 2.0 * fy[μi, σi]) 
            v∂μ = 0.5 * (fy[μi+1, σi] - fy[μi-1, σi]) 
        elseif μi == 1
            v∂μ = - 0.5 * (3 * fy[μi, σi] - 4 * fy[μi+1, σi] + fy[μi+2, σi])
        elseif μi == μn
            v∂μ = 0.5 * (3 * fy[μi, σi] - 4 * fy[μi-1, σi] + fy[μi-2, σi])
        end
        ydot[ij] = (byp.ρ * byp.θ * max(fy[μi, σi], 0.0)^(1-1/byp.θ)
                    + (- byp.ρ * byp.θ + (1-byp.γ) * byp.μs[μi] - 0.5 * (1-byp.γ) * byp.γ * byp.νD^2 * byp.σs[σi]) * fy[μi, σi] 
                    + ∂μ*v∂μ 
                    + ∂σ*v∂σ 
                    + ∂2μ*v∂2μ 
                    + ∂2σ*v∂2σ 
                    )
    end
    return ydot
end

function F(byp::BansalYaronProblem, y::Vector)
    ydot = deepcopy(y)
    F!(byp, y, ydot)
end

function solve(::Type{Val{:ode45}}, byp::BansalYaronProblem; kwargs...)
    oldV = deepcopy(byp.V)
    newV = oldV
    distance = Inf
    i = 0
    while i < 10
        i += 1
        @show i, distance
        xout, yout = ODE.ode45((t, y) -> (ydot = deepcopy(y) ; F!(byp, y, ydot) ; return y), oldV, [0.0;100.0]; kwargs...)
        newV = yout[end]
        distance = chebyshev(newV, oldV)
        if distance < 1e-5
            break
        else
            oldV = newV
        end
    end
    return newV
end


function solve(::Type{Val{:adiff}}, byp::BansalYaronProblem; kwargs...)
    oldV = deepcopy(byp.V)
    df = DifferentiableMultivariateFunction((y, ydot) -> F!(byp, y, ydot))
    out = nlsolve(df, oldV, method = :trust_region, show_trace = true)
    return out.zero
end


dμ=Interval(min(byp.μs),max(byp.μs))
dσ=Interval(min(byp.σs),max(byp.σs))
d = dμ * dσ
Dμ=Derivative(d,1)
Dσ=Derivative(d,2)

D2μ=Derivative(d,2)
D2σt=Derivative(d,2)

;dt=Interval(0.0,.1)
d=dx*dt

V=Fun(x->x^2,dx)

Dt=Derivative(d,2);Dx=Derivative(d,1)

ϵ=1.
