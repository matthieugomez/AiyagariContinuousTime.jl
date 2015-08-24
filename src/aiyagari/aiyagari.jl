# requires Julia 0.4
##############################################################################
##
## Type AiygariProblem
##
##############################################################################
typealias sMatrix{Tv, Ti} Base.SparseMatrix.SparseMatrixCSC{Tv, Ti}

type AiyagariProblem
    π::Float64               # productivity
    γ::Float64               # CRRA utility with parameter gamma
    α::Float64               # Production function F = K^α * L^(1-α) 
    δ::Float64               # Capital depreciation
    ρ::Float64               # discount rate
    
    zn::Int                  # number of productivity points 
    zmin::Float64            # min productivity
    zmax::Float64            # max productivity
    z::Vector{Float64}       # productivity vector
  
    an::Int                  # number of asset points 
    amin::Float64            # min asset
    amax::Float64            # max asset
    a::Vector{Float64}       # wealth vector
  
    invΔ::Float64            # 1/δ in HznB algorithm

    C::sMatrix{Float64, Int} # C matrix (transposed)
   
    b::Vector{Float64}
   
    K::Float64               # initial aggregate capital. 
    V::Vector{Float64}       # initial value function
  
    gg::Vector{Float64}      # distribution
    newV::Vector{Float64}    # old value function


    # storage arrays
    A::sMatrix{Float64,Int}   # An matrix (transposed)
    B::sMatrix{Float64,Int}   # B matrix (transposed)
    r::Float64
    w::Float64
    ra::Vector{Float64}
    wz::Vector{Float64}
    ∂Vf::Matrix{Float64}
    ∂Vb::Matrix{Float64}
    ∂V0::Matrix{Float64}
    sf::Matrix{Float64}
    sb::Matrix{Float64}
    u::Vector{Float64}
end

# construt an instance of AiyagariProblem
function AiyagariProblem(π::Float64 = 1.0, 
                         K::Float64 = 3.8; 
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
                         invΔ::Float64 = 1e-3
                         )

    a = collect(linspace(amin,amax,an))  
    z = collect(linspace(zmin,zmax,zn))  
    invdz = (zn-1) / (zmax-zmin)

    # create a sparse matrix of the right form
    A = spdiagm(
        (ones(an*zn), ones(an*zn-1), ones(an*zn-1), ones((zn-1)*an), ones((zn-1)*an)),
        (0, 1, -1, an, -an)
    )
    B = deepcopy(A)
    C = deepcopy(A)

    # fill up matrix C
    Cvals = nonzeros(C) 
    Crows = rowvals(C) 
    fill!(Cvals, zero(Float64))
    ij = zero(Int)
    @inbounds for zi in 1:zn, ai in 1:an
        ij += 1
        krange = nzrange(C, ij)
        rows = Crows[krange]

        current =  0.5 * σ2 * invdz^2
        if zi > 1
            index = searchsortedfirst(rows, ij - an)
        else
            index = searchsortedfirst(rows, ij)
        end
        Cvals[krange[index]] += current

        current = θ * (zmean - z[zi]) * invdz + 0.5 * σ2 * invdz^2
        if zi < zn
            index = searchsortedfirst(rows, ij + an)
        else
            index = searchsortedfirst(rows, ij)
        end
        Cvals[krange[index]] += current

        current = - σ2 * invdz^2 - θ * (zmean - z[zi]) * invdz
        index = searchsortedfirst(rows, ij)
        Cvals[krange[index]] += current
    end

    # b such that Ag = b in plank
    b = fill(zero(Float64), an*zn)
    i_fix = 1
    b[i_fix] = .1

    # initial value function
    r = π * α * K^(α-1) - δ 
    w = π * (1-α) * K^α   
    V = Array(Float64, an*zn)
    index = zero(Int)
    for zi in 1:zn, ai in 1:an
        index += 1
        V[index] = (w * z[zi] + r * a[ai])^(1-γ)
    end
    scale!(V, 1/((1-γ)*ρ))
    newV = deepcopy(V)

    # storage arrays
    gg = Array(Float64, an*zn)
    ∂Vf = Array(Float64, an, zn)
    ∂Vb = Array(Float64, an, zn)
    ∂V0 = Array(Float64, an, zn)
    sf = Array(Float64, an, zn)
    sb = Array(Float64, an, zn)
    u = Array(Float64, an*zn)
    wz = similar(z)
    ra = similar(a)

    AiyagariProblem(π, γ, α, δ, ρ, zn, zmin, zmax, z, an, amin, amax, a, invΔ, C, b, K, V, gg, newV, A, B, r, w, ra, wz, ∂Vf, ∂Vb, ∂V0, sf, sb, u)
end





##############################################################################
##
## Utilities functions used in solving both  the stationary and dynamic problem
##
##############################################################################


function update_value!(ap::AiyagariProblem)

    amin = ap.amin ; amax = ap.amax ; γ = ap.γ ; an = ap.an ; zn = ap.zn ; 
    invΔ = ap.invΔ ;  ρ = ap.ρ ; C = ap.C ; wz = ap.wz ; ra = ap.ra ; A = ap.A ;
    z = ap.z ; a = ap.a ; w = ap.w ; r = ap.r ;  ∂Vf = ap.∂Vf; ∂Vb = ap.∂Vb ; 
    ∂V0 = ap.∂V0 ; sf = ap.sf ; sb = ap.sb ; u = ap.u ; B = ap.B ; V = ap.V ; 

    # precompute some quantities
    invda::Float64 = (an-1)/(amax-amin)
    invγ::Float64 = 1/γ
    anzn::Int = an*zn
    inv1γ::Float64 = 1/(1-γ)
    ramax::Float64 = r * amax
    ramin::Float64 = r * amin

    # update ∂Vf (forward)
    index = zero(Int)
    @inbounds for zi in 1:zn
        for ai in 1:(an-1)
            index += 1
            ∂Vf[ai, zi] = (V[index+1] - V[index]) * invda
        end
        # state constraint boundary condition
        # will never be used, but just in case
        index += 1
        ∂Vf[an, zi] = (wz[zi] + ramax)^(-γ)
    end

    # update ∂Vb (backward)
    index = zero(Int)
    @inbounds for zi in 1:zn
        # state constraint boundary condition
        index += 1
        ∂Vb[1, zi] = (wz[zi] + ramin)^(-γ)  
        for ai in 2:an
            index += 1
            ∂Vb[ai, zi] = (V[index] - V[index-1]) * invda
        end
    end

    # update ∂V0,  derivative of value function at steady state
    @inbounds for zi in 1:zn
        wzi = wz[zi]
        @simd for ai in 1:an
            ∂V0[ai, zi] = (wzi + ra[ai])^(-γ)
        end
    end

    # update sf
    @inbounds for zi in 1:zn
        @simd for ai in 1:an
            sf[ai, zi] = wz[zi] + ra[ai] - ∂Vf[ai, zi]^(-invγ)
        end
    end

    # update sb
    @inbounds for zi in 1:zn
        @simd for ai in 1:an
            sb[ai, zi] = wz[zi] + ra[ai] - ∂Vb[ai, zi]^(-invγ)
        end
    end

    # update u
    ij = zero(Int)
    @inbounds for zi in 1:zn, ai in 1:an
        ij += 1
        value = sf[ai, zi] > 0  ? ∂Vf[ai, zi] : sb[ai, zi] < 0 ? ∂Vb[ai, zi] : ∂V0[ai, zi]
        u[ij] = inv1γ * value^(1-invγ) + invΔ * V[ij]
    end

    # update A
    Cvals = nonzeros(C)
    Avals = nonzeros(A)
    copy!(Avals, Cvals)
    Arows = rowvals(A) 
    ij = zero(Int)
    @inbounds for zi in 1:zn, ai in 1:an
        ij += 1
        krange = nzrange(A, ij)
        rows = Arows[krange]

        current = -min(sb[ai, zi], 0.0) * invda 
        if ai > 1
            index = searchsortedfirst(rows, ij - 1)
        else
            index = searchsortedfirst(rows, ij)
        end
        Avals[krange[index]] += current

        current = max(sf[ai, zi], 0.0) * invda 
        if ai < an
            index = searchsortedfirst(rows, ij + 1)
        else
            index = searchsortedfirst(rows, ij)
        end
        Avals[krange[index]] += current

        current = min(sb[ai, zi], 0.0) * invda - max(sf[ai, zi], 0.0) * invda 
        index = searchsortedfirst(rows, ij)
        Avals[krange[index]] += current
    end

    # update B = diag(Δ + ρ) - A
    index = zero(Int)
    Brows = rowvals(B)
    Avals = nonzeros(A) 
    Bvals = nonzeros(B) 
    ij = zero(Int)
    @inbounds for zi in 1:zn, ai in 1:an
        ij += 1
        for k in nzrange(B, ij)
            # loop over elements in the column ij
            row = Brows[k]
            Bvals[k] = - Avals[k]
            if row == ij
                Bvals[k] += invΔ + ρ
            end
        end
    end

    # solve for new V
    ap.newV = B' \ u
end

function chebyshev{N}(x1::Array{Float64, N}, x2::Array{Float64, N})
    out = zero(Float64)
    @inbounds @simd for i in 1:length(x1)
        current = abs(x1[i] - x2[i])
        if current > out
            out = current
        end
    end
    return out
end


function sum_capital(a::Vector{Float64}, gg::Vector{Float64})
    numerator = zero(Float64)
    denominator = zero(Float64)
    an = length(a)
    zn = div(length(gg), length(a))
    @inbounds for zi in 1:zn
        start = (zi-1)*an
        @simd for ai in 1:an
            current = gg[start + ai]
            numerator += current * a[ai] 
            denominator += current
        end
    end
    return numerator / denominator
end


##############################################################################
##
## Solve static problem
##
##############################################################################

# solve hmj
function solve_hmj!(ap::AiyagariProblem ;  
                    maxit::Int = 100, 
                    crit::Float64 = 1e-6, 
                    verbose::Bool = true)
    # update price vectors
    ap.r = ap.π * ap.α * ap.K^(ap.α-1) - ap.δ    
    ap.w = ap.π * (1-ap.α) * ap.K^ap.α   
    broadcast!(*, ap.wz, ap.z, ap.w)
    broadcast!(*, ap.ra, ap.a, ap.r)
    for iter in 1:maxit
        # update newV
        update_value!(ap)
        # check convergence
        distance = chebyshev(ap.newV, ap.V)
        if distance < crit
            if verbose
                println("hmj solved : $(iter) iterations")
            end
            break
        else
            # update V using newV
            (ap.newV, ap.V) = (ap.V, ap.newV)
        end
    end
end

# solve fokker-planck
function solve_fp!(ap::AiyagariProblem)
    A = ap.A 
    # change first row of matrix A
    Arows = rowvals(A)
    Avals = nonzeros(A)
    @inbounds for ij in 1:size(A, 2)
        for k in nzrange(A, ij)
            if Arows[k] == 1
                Avals[k] = zero(Float64)
            end
        end
    end
    Avals[1] = one(Float64)
    
    # solve system Ag = b
    ap.gg = A \ ap.b
end


# solve
function solve!(ap::AiyagariProblem; 
                maxit::Int  = 100,      # maximum number of iterations in the HznB loop
                maxitK::Int = 100,      # maximum number of iterations in the K loop
                crit::Float64 = 1e-6,   # criterion HznB loop
                critK::Float64 = 1e-5,  # criterion K loop
                relax::Float64 = 0.99,  # relaxation parameter 
                verbose = true)
    for iter in 1:maxitK
        if verbose
            println("Main loop iteration ",iter)       
        end

        # solve hamilton-jacobi-bellman equation 
        solve_hmj!(ap, maxit = maxit, crit = crit, verbose = verbose)

        # solve fokker-planck equation
        solve_fp!(ap)

        # Update aggregate capital
        newK = sum_capital(ap.a, ap.gg)
        if verbose
            @show newK
        end
        if abs(ap.K - newK) < critK
            break
        else
            # relaxation algorithm (to ensure convergence)
            ap.K = relax * ap.K + (1 - relax) * newK  
        end
    end
end

##############################################################################
##
## Type DynamicAiygariProblem
##
##############################################################################

type DynamicAiyagariProblem
    ap::AiyagariProblem 
    π::Vector{Float64}               # productivity across time
    K::Vector{Float64}               # capital across time
    Vend::Vector{Float64}            # value function at the end
    gg::Vector{Vector{Float64}}      # distribution. gg[1] is initial condition.
    A::Vector{sMatrix{Float64, Int}} # A matrix across time
    N::Int                           # Number of periods
    dt::Float64
end

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
## Solve Dynamic Problem
##
##############################################################################

# solve hmj
function solve_hmj!(dap::DynamicAiyagariProblem)
    dap.ap.V = dap.Vend
    invdt = 1 / dap.dt
    for n in dap.N:-1:1
        # create subproblem
        ap = AiyagariProblem(dap.π[n], dap.ap.γ, dap.ap.α, dap.ap.δ, dap.ap.ρ, dap.ap.zn, dap.ap.zmin, dap.ap.zmax, dap.ap.z, dap.ap.an, dap.ap.amin, dap.ap.amax, dap.ap.a, invdt, dap.ap.C, dap.ap.b, dap.K[n], dap.ap.V, dap.gg[n], dap.ap.newV, dap.A[n], dap.ap.B, dap.ap.r, dap.ap.w, dap.ap.ra, dap.ap.wz, dap.ap.∂Vf, dap.ap.∂Vb, dap.ap.∂V0, dap.ap.sf, dap.ap.sb, dap.ap.u)

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
        # convert A to I - A
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
        println("Main loop iteration ",iter)
        solve_hmj!(dap)
        solve_fp!(dap)
        # Update aggregate capital
        for n in 1:N
            newK[n] = sum_capital(dap.ap.a, dap.gg[n])
        end
        distance = chebyshev(dap.K, newK)
        @show distance
        if chebyshev(dap.K, newK) < critK
            break
        else
            # relaxation algorithm (to ensure convergence)
            dap.K = relax * dap.K .+ (1 - relax) * newK  
        end
    end
end

