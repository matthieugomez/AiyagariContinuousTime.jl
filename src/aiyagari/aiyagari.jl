#  Achdou, Han, Lasry, Lions, Moll (2015)
#  Heterogeneous Agent Models in Continuous Time


##############################################################################
##
## Type
##
##############################################################################

type AiyagariProblem
    π::Float64               # productivity
    α::Float64               # Production function F = K^α * L^(1-α) 
    δ::Float64               # Capital depreciation
    γ::Float64               # CRRA utility with parameter gamma
    ρ::Float64               # discount rate
    σ2::Float64
    θ::Float64
end

function AiyagariProblem(; π::Float64 = 1.0,   α::Float64 = 0.35, δ::Float64 = 0.1, γ::Float64 = 2.0, ρ::Float64 = 0.05, σ2::Float64 = (0.10)^2,   θ::Float64 = 0.3)
    AiyagariProblem(π, α, δ, γ , ρ, σ2, θ)
end

abstract AiyagariMethod
##############################################################################
##
## Solve equilibrium by iterating on K
##
##############################################################################

# compute aggregate capital
function sum_capital(a::Vector{Float64}, gg::Vector{Float64})
    numerator = zero(Float64)
    denominator = zero(Float64)
    an = length(a)
    zn = div(length(gg), length(a))
    ij = zero(Int)
    @inbounds for zi in 1:zn
        @simd for ai in 1:an
            ij += 1
            current = gg[ij]
            numerator += current * a[ai] 
            denominator += current
        end
    end
    return numerator / denominator
end

# solve
function solve!(afd::AiyagariMethod; 
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
        solve_hjb!(afd, maxit = maxit, crit = crit, verbose = verbose)
        v = afd.V[(length(afd.a) * 10 + 1):(length(afd.a) * 11)]
        @show v
        # solve fokker-planck equation
        solve_fp!(afd)
        g = afd.gg[1:length(afd.a)]
        @show g
        # Update aggregate capital
        newK = sum_capital(afd.a, afd.gg)
        if verbose
            @show newK
        end
        if abs(afd.K - newK) < critK
            break
        else
            # relaxation algorithm (to ensure convergence)
            afd.K = relax * afd.K + (1 - relax) * newK  
        end
    end
end



##############################################################################
##
## Constructor
##
##############################################################################

type AiyagariFD <: AiyagariMethod
    ap::AiyagariProblem 

    # discretization
    z::Vector{Float64}       # productivity vector
    a::Vector{Float64}       # wealth vector
    invΔ::Float64            # 1/δ in HznB algorithm

    # solution
    K::Float64               # aggregate capital 
    r::Float64
    w::Float64
    V::Vector{Float64}       # ivalue function
    gg::Vector{Float64}      # distribution
    A::Base.SparseMatrix.SparseMatrixCSC{Float64,Int}   

    # storage 
    b::Vector{Float64}
    C::Base.SparseMatrix.SparseMatrixCSC{Float64, Int} 
    newV::Vector{Float64}    # old value function
    B::Base.SparseMatrix.SparseMatrixCSC{Float64,Int}   
    u::Vector{Float64}
end



function AiyagariFD(ap::AiyagariProblem,  
                     zn::Int = 30, 
                     zmin::Float64 = 0.5, 
                     zmax::Float64 = 1.5, 
                     an::Int = 30, 
                     amin::Float64 = -1.0, 
                     amax::Float64 = 30.0, 
                     invΔ::Float64 = 1e-3,
                     K::Float64 = 3.8
                     )
    π = ap.π ; α = ap.α ; δ = ap.δ ; γ = ap.γ ; ρ = ap.ρ ;  σ2 = ap.σ2 ; θ = ap.θ

    a = collect(linspace(amin,amax,an))  
    z = collect(linspace(zmin,zmax,zn))  
    invdz = 1/ (z[2]-z[1])

    # create a sparse matrix of the right form filled with zeros
    A = spdiagm(
        (ones(an*zn), ones(an*zn-1), ones(an*zn-1), ones((zn-1)*an), ones((zn-1)*an)),
        (0, 1, -1, an, -an)
    )
    fill!(nonzeros(A), zero(Float64))
    B = deepcopy(A)
    C = deepcopy(A)

    # fill up matrix C
    Cvals = nonzeros(C) 
    Crows = rowvals(C) 
    ij = zero(Int)
    @inbounds for zi in 1:zn, ai in 1:an
        ij += 1
        if zi > 1
            current = 0.5 * σ2 * invdz^2
            C[ij - an, ij] += current
            C[ij, ij] -= current
        end    
        if zi < zn
            current = θ * (1.0 - z[zi]) * invdz + 0.5 * σ2 * invdz^2
            C[ij + an, ij] += current
            C[ij, ij] -= current
        end
    end

    # b such that Ag = b in plank
    b = fill(zero(Float64), an*zn)
    i_fix = 1
    b[i_fix] = .1

    # initial value function
    r = π * α * K^(α-1) - δ 
    w = π * (1-α) * K^α   
    V = Array(Float64, an*zn)
    ij = zero(Int)
    for zi in 1:zn, ai in 1:an
        ij += 1
        V[ij] = (w * z[zi] + r * a[ai])^(1-γ)/((1-γ)*ρ)
    end
    newV = deepcopy(V)

    # storage arrays
    gg = Array(Float64, an*zn)
    u = Array(Float64, an*zn)

    AiyagariFD(ap, z, a, invΔ, K, r, w, V, gg, A, b, C, newV, B, u)
end


##############################################################################
##
## Solve hamilton-jacobi-bellman
##
##############################################################################

function update_value!(afd::AiyagariFD)
    π = afd.ap.π ; α = afd.ap.α ; δ = afd.ap.δ ; γ = afd.ap.γ ; ρ = afd.ap.ρ ;  σ2 = afd.ap.σ2 ; θ = afd.ap.θ
    invΔ = afd.invΔ  ; C = afd.C ; A = afd.A ; z = afd.z ; a = afd.a ; w = afd.w ; r = afd.r   ; u = afd.u ; B = afd.B ; V = afd.V ; 

    # precompute some quantities
    an = length(a)
    zn = length(z)
    invda = 1/(a[2]-a[1])
    invγ = 1/γ
    inv1γ = 1.0/(1.0-γ)

    # set A = C
    Cvals = nonzeros(C)
    Avals = nonzeros(A)
    copy!(Avals, Cvals)

    # update A and u
    ij = zero(Int)
    @inbounds for zi in 1:zn, ai in 1:an
        ij += 1
        krange = nzrange(A, ij)
        rows = rowvals(A)[krange]
        if ai < an
            ∂V = (V[ij+1] - V[ij]) * invda
        else
            # state constraint boundary condition
            ∂V = (z[zi] * afd.w + a[end] * afd.r)^(-γ)
        end
        saving = z[zi] * afd.w + a[ai] * afd.r - ∂V^(-invγ)
        if saving > 0 # case of positive drift
            if ai < an
                current = saving * invda 
                A[ij, ij] -= current
                A[ij + 1, ij] += current
            end
        else
            if ai > 1
                ∂V = (V[ij] - V[ij-1]) * invda
            else
                # state constraint boundary condition
                ∂V =  (z[zi] * afd.w + a[1] * afd.r)^(-γ) 
            end
            saving = z[zi] * afd.w + a[ai] * afd.r - ∂V^(-invγ)
            if saving < 0 # case of negative drift
                if ai > 1
                    current = saving * invda 
                    A[ij - 1, ij] -= current
                    A[ij, ij] += current
                end
            else
                ∂V = (z[zi] * afd.w + a[ai] * afd.r)^(-γ)
            end
        end

        # update u
        u[ij] = inv1γ * ∂V^(1-invγ) + invΔ * V[ij]
    end
    
    # set B = diag(Δ + ρ) - A
    Bvals = nonzeros(B) 
    Brows = rowvals(B)
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
    afd.newV = B' \ u
end

function solve_hjb!(afd::AiyagariFD ;  
                    maxit::Int = 100, 
                    crit::Float64 = 1e-6, 
                    verbose::Bool = true)
    # update price vectors
    afd.r = afd.ap.π * afd.ap.α * afd.K^(afd.ap.α-1) - afd.ap.δ    
    afd.w = afd.ap.π * (1-afd.ap.α) * afd.K^afd.ap.α   
    for iter in 1:maxit
        # update newV
        update_value!(afd)
        # check convergence
        distance = chebyshev(vec(afd.newV), vec(afd.V))
        if distance < crit
            if verbose
                println("hjb solved : $(iter) iterations")
            end
            break
        else
            # update V using newV
            (afd.newV, afd.V) = (afd.V, afd.newV)
        end
    end
end


##############################################################################
##
## Solve fokker-planck equation
##
##############################################################################

# solve fokker-planck
function solve_fp!(afd::AiyagariFD)
    A = afd.A
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
    afd.gg  = A \ afd.b
end

