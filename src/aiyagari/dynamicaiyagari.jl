#  Achdou, Han, Lasry, Lions, Moll (2015)
#  Heterogeneous Agent Models in Continuous Time

##############################################################################
##
## Type 
##
##############################################################################

type DynamicAiyagariProblem
    ap::AiyagariProblem 
    π::Vector{Float64}               # productivity across time
    K::Vector{Float64}               # capital across time
    Vend::Vector{Float64}            # value function at the end
    gg::Vector{Vector{Float64}}      # distribution. gg[1] is initial condition.
    A::Vector{Base.SparseMatrix.SparseMatrixCSC{Float64, Int}} # A matrix across time
    N::Int                           # Number of periods
    dt::Float64
end

##############################################################################
##
## Constructor
##
##############################################################################

# Construct an instance of the tyoe if initial values for K is a number
function DynamicAiyagariProblem(
    π::Vector{Float64}, 
    K::Float64; 
    dt::Float64 = 0.5, 
    γ::Float64 = 2.0, 
    α::Float64 = 0.35, 
    δ::Float64 = 0.1, 
    ρ::Float64 = 0.05, 
    σ2::Float64 = (0.10)^2,  
    θ::Float64 = 0.3,
    zmean::Float64 = 1.0,      
    zn::Int = 40, 
    zmin::Float64 = 0.5, 
    zmax::Float64 = 1.5, 
    amin::Float64 = -1.0, 
    amax::Float64 = 30.0, 
    an::Int = 100, 
    invΔ::Float64 = 1e-3, 
    relax::Float64 = 0.99)

    # initialize path K by solving stationary model at every π
    ap = AiyagariProblem(π[1], K, γ = γ, α = α, δ = δ, ρ = ρ, σ2 = σ2, θ = θ, zmean = zmean, zn = zn, zmin = zmin, zmax = zmax, amin = amin, amax = amax, an = an, invΔ = invΔ)
    N = length(π)
    Ks = Array(Float64, N)
    N0 = floor(Int, N/10)

    for n in 1:N0
        ap.π = π[n]
        if n > 1
            ap.K = Ks[n-1]
        end
        solve!(ap, relax = relax, verbose = false)
        Ks[n] = ap.K
    end
    for n in (N0+1):N
        Ks[n] = ap.K
    end

    DynamicAiyagariProblem(π, Ks, dt = dt, γ = γ, α = α, δ = δ, ρ = ρ, σ2 = σ2, θ = θ, zmean = zmean, zn = zn, zmin = zmin, zmax = zmax, amin = amin, amax = amax, an = an, invΔ = invΔ, relax = relax)
end

# Construct an instance of the tyoe if initial values for K is a vector
# K[t] should value of stationary equilibrium π[t]
function  DynamicAiyagariProblem(
    π::Vector{Float64}, 
    K::Vector{Float64}; 
    dt::Float64 = 0.5, 
    γ::Float64 = 2.0, 
    α::Float64 = 0.35, 
    δ::Float64 = 0.1, 
    ρ::Float64 = 0.05, 
    σ2::Float64 = (0.10)^2,  
    θ::Float64 = 0.3,
    zmean::Float64 = 1.0,    
    zn::Int = 40, 
    zmin::Float64 = 0.5, 
    zmax::Float64 = 1.5, 
    amin::Float64 = -1.0, 
    amax::Float64 = 30.0, 
    an::Int = 100,
    invΔ::Float64 = 1e-3, 
    relax::Float64 = 0.99)

    ap = AiyagariProblem(π[end], K[end], γ = γ, α = α, δ = δ, ρ = ρ, σ2 = σ2, θ = θ, zmean = zmean, zn = zn, zmin = zmin, zmax = zmax, amin = amin, amax = amax, an = an, invΔ = invΔ)
    solve!(ap, relax = relax, verbose = false)
    N = length(π)

    # initialize storage
    A = Array(Base.SparseMatrix.SparseMatrixCSC{Float64,Int64}, N)
    gg = Array(Vector{Float64}, N)
    for n in 1:N
        A[n] = similar(ap.A)
        gg[n] = similar(ap.gg)
    end

    # boundary condition 1: value function at n = N
    Vend = deepcopy(ap.V)
    ## boundary condition 1: distribution at n = 1
    gg[1] = deepcopy(ap.gg)


    DynamicAiyagariProblem(ap, π, K, Vend, gg, A, N, dt)
end


##############################################################################
##
## Solve 
##
##############################################################################

# solve hjb
function solve_hjb!(dap::DynamicAiyagariProblem)
    dap.ap.V = dap.Vend
    invdt = 1 / dap.dt
    for n in dap.N:-1:1
        # create static problem at date n
        ap = AiyagariProblem(dap.π[n], dap.ap.γ, dap.ap.α, dap.ap.δ, dap.ap.ρ, dap.ap.z, dap.ap.a, invdt, dap.ap.C, dap.ap.b, dap.K[n], dap.ap.V, dap.gg[n], dap.ap.newV, dap.A[n], dap.ap.B, dap.ap.u, dap.ap.r, dap.ap.w, dap.ap.ra, dap.ap.wz)

        # update prices
        ap.r = ap.π * ap.α * ap.K^(ap.α-1) - ap.δ    
        ap.w = ap.π * (1-ap.α) * ap.K^ap.α   
        broadcast!(*, ap.wz, ap.z, ap.w)
        broadcast!(*, ap.ra, ap.a, ap.r)

        # update value
        update_value!(ap)

        (dap.ap.newV, dap.ap.V) = (dap.ap.V, dap.ap.newV)      
    end
end

# solve fokker-planck
function solve_fp!(dap::DynamicAiyagariProblem)
    for n in 1:(dap.N-1)
        A = dap.A[n]
        # convert A to I - dt * A
        Avals = nonzeros(A) 
        Arows = rowvals(A)
        @inbounds for ij in 1:size(A, 2)
            for k in nzrange(A, ij)
                Avals[k] = - dap.dt * Avals[k]
                if Arows[k] == ij
                    Avals[k] += one(Float64)
                end
            end
        end
        # obtain next period distribution
        dap.gg[n+1] = A \ dap.gg[n]
    end
end

function solve!(dap::DynamicAiyagariProblem ;
                maxitK::Int = 100,            # maximum number of iterations in the K loop
                critK::Float64 = 1e-5,        # criterion K loop
                relax::Float64 = 0.1,         # relaxation parameter 
                )
    newK = similar(dap.K)
    for iter in 1:maxitK
        println("Main loop iteration ", iter)
        solve_hjb!(dap)
        solve_fp!(dap)
        # Update aggregate capital
        for n in 1:dap.N
            newK[n] = sum_capital(dap.ap.a, dap.gg[n])
        end
        distance = chebyshev(dap.K, newK)
        @show distance
        if chebyshev(dap.K, newK) < critK
            break
        else
            # relaxation algorithm 
            dap.K = relax * dap.K .+ (1 - relax) * newK  
        end
    end
end
