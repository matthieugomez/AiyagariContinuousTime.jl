#  Achdou, Han, Lasry, Lions, Moll (2015)
#  Heterogeneous Agent Models in Continuous Time


##############################################################################
##
## Type
##
##############################################################################

type AiyagariProblem
    π::Float64               # productivity
    γ::Float64               # CRRA utility with parameter gamma
    α::Float64               # Production function F = K^α * L^(1-α) 
    δ::Float64               # Capital depreciation
    ρ::Float64               # discount rate
    
    z::Vector{Float64}       # productivity vector
    a::Vector{Float64}       # wealth vector
  
    invΔ::Float64            # 1/δ in HznB algorithm

    C::Base.SparseMatrix.SparseMatrixCSC{Float64, Int} 
   
    b::Vector{Float64}
   
    K::Float64               # initial aggregate capital. 
    V::Vector{Float64}       # initial value function
  
    gg::Vector{Float64}      # distribution
    newV::Vector{Float64}    # old value function


    # storage 
    A::Base.SparseMatrix.SparseMatrixCSC{Float64,Int}   
    B::Base.SparseMatrix.SparseMatrixCSC{Float64,Int}   
    u::Vector{Float64}
    r::Float64
    w::Float64
    ra::Vector{Float64}
    wz::Vector{Float64}
end

##############################################################################
##
## Constructor
##
##############################################################################


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
        krange = nzrange(C, ij)
        rows = Crows[krange]

        if zi > 1
            current =  0.5 * σ2 * invdz^2
            index = searchsortedfirst(rows, ij - an)
            Cvals[krange[index]] += current
            index = searchsortedfirst(rows, ij)
            Cvals[krange[index]] -= current
        end    
        
        if zi < zn
            current = θ * (zmean - z[zi]) * invdz + 0.5 * σ2 * invdz^2
            index = searchsortedfirst(rows, ij + an)
            Cvals[krange[index]] += current
            index = searchsortedfirst(rows, ij)
            Cvals[krange[index]] -= current
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
    wz = similar(z)
    ra = similar(a)

    AiyagariProblem(π, γ, α, δ, ρ, z, a, invΔ, C, b, K, V, gg, newV, A, B, u, r, w, ra, wz)
end


##############################################################################
##
## Solve hamilton-jacobi-bellman
##
##############################################################################

function update_value!(ap::AiyagariProblem)
    γ = ap.γ  ;  a = ap.a ; z = ap.z
    invΔ = ap.invΔ ;  ρ = ap.ρ ; C = ap.C ; wz = ap.wz ; ra = ap.ra ; A = ap.A ;
    z = ap.z ; a = ap.a ; w = ap.w ; r = ap.r   ; u = ap.u ; B = ap.B ; V = ap.V ; 

    # precompute some quantities
    an = length(a)
    zn = length(z)
    invda = (an-1)/(a[end]-a[1])
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
            ∂V = (wz[zi] + ra[end])^(-γ)
        end
        saving = wz[zi] + ra[ai] - ∂V^(-invγ)
        if saving > 0 # case of positive drift
            if ai < an
                current = saving * invda 
                index = searchsortedfirst(rows, ij)
                Avals[krange[index]] -= current
                index = searchsortedfirst(rows, ij + 1)
                Avals[krange[index]] += current
            end
        else
            if ai > 1
                ∂V = (V[ij] - V[ij-1]) * invda
            else
                # state constraint boundary condition
                ∂V =  (wz[zi] + ra[1])^(-γ) 
            end
            saving = wz[zi] + ra[ai] - ∂V^(-invγ)
            if saving < 0 # case of negative drift
                if ai > 1
                    current = saving * invda 
                    index = searchsortedfirst(rows, ij - 1)
                    Avals[krange[index]] -= current
                    index = searchsortedfirst(rows, ij)
                    Avals[krange[index]] += current
                end
            else
                ∂V = (wz[zi] + ra[ai])^(-γ)
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
    ap.newV = B' \ u
end

function solve_hjb!(ap::AiyagariProblem ;  
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
        distance = chebyshev(vec(ap.newV), vec(ap.V))
        if distance < crit
            if verbose
                println("hjb solved : $(iter) iterations")
            end
            break
        else
            # update V using newV
            (ap.newV, ap.V) = (ap.V, ap.newV)
        end
    end
end


##############################################################################
##
## Solve fokker-planck equation
##
##############################################################################

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
        solve_hjb!(ap, maxit = maxit, crit = crit, verbose = verbose)

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