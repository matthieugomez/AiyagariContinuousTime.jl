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
    σ2::Float64
    θ::Float64
    zmean::Float64
    z::Vector{Float64}       # productivity vector
    a::Vector{Float64}       # wealth vector
    invΔ::Float64            # 1/δ in HznB algorithm
end

function AiyagariProblem(;
                         π::Float64 = 0.0, 
                         γ::Float64 = 2.0, 
                         α::Float64 = 0.35, 
                         δ::Float64 = 0.1, 
                         ρ::Float64 = 0.05, 
                         σ2::Float64 = (0.10)^2,  
                         θ::Float64 = 0.3,
                         zmean::Float64 = 1.0,      
                         zn::Int = 20, 
                         zmin::Float64 = 0.5, 
                         zmax::Float64 = 1.5, 
                         amin::Float64 = -1.0, 
                         amax::Float64 = 50.0, 
                         an::Int = 50, 
                         invΔ::Float64 = 1e-3
                         )

    a = collect(linspace(amin,amax,an))  
    z = collect(linspace(zmin,zmax,zn))  
    AiyagariProblem(π, γ, α, δ, ρ, σ2, θ, zmean, z, a, invΔ)
end



##############################################################################
##
## Aiyagari Solution
##
##############################################################################
type AiyagariSolution{T}
    V::Vector{T}       # initial value function
    g::Vector{T}      # distribution
    K::T               # initial aggregate capital. 
    r::T
    w::T
end

function AiyagariSolution{T}(ap::AiyagariProblem; K::T = 3.8)
    # initial value function
    π = ap.π ; α = ap.α ; δ = ap.δ ; z = ap.z ; a = ap.a ; γ = ap.γ ; ρ = ap.ρ
    an = length(a)
    zn = length(z)
    r = exp(π) * α * K^(α-1) - δ 
    w = exp(π) * (1-α) * K^α   
    V = Array(T, an*zn)
    ij = zero(Int)
    for zi in 1:zn, ai in 1:an
        ij += 1
        V[ij] = (w * z[zi] + r * a[ai])^(1-γ)/((1-γ)*ρ)
    end

    # storage arrays
    g = Array(T, an*zn)
    AiyagariSolution(V, g, K, r, w)
end

##############################################################################
##
## Aiyagari Storage
##
##############################################################################

type AiyagariArrays{T}
    C::Base.SparseMatrix.SparseMatrixCSC{T, Int} 
    B::Base.SparseMatrix.SparseMatrixCSC{T, Int}   
    A::Base.SparseMatrix.SparseMatrixCSC{T, Int}   
    b::Vector{T}
    u::Vector{T}
end

function AiyagariArrays{T}(ap::AiyagariProblem, as::AiyagariSolution{T})
    π = ap.π ; α = ap.α ; δ = ap.δ ; z = ap.z ; a = ap.a ; γ = ap.γ ; ρ = ap.ρ; σ2 = ap.σ2 ; θ = ap.θ ; zmean = ap.zmean
    invdz = 1/(z[2] - z[1])
    invda = 1/(a[2] - a[1])
    an = length(a)
    zn = length(z)

    # create 3 sparse matrices of the right form filled with zeros
    A = spdiagm(
        (ones(T, an*zn), ones(T, an*zn-1), ones(T, an*zn-1), ones(T, (zn-1)*an), ones(T, (zn-1)*an)),
        (0, 1, -1, an, -an)
    )
    fill!(nonzeros(A), zero(T))
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
    b = fill(zero(T), an*zn)
    b[1] = .1

    u = similar(b)

    AiyagariArrays(C, B, A, b, u)
end


##############################################################################
##
## Solve hamilton-jacobi-bellman
##
##############################################################################

    # update A and u
function  update_Au!(ap::AiyagariProblem, aa::AiyagariArrays, as::AiyagariSolution)
    γ = ap.γ  ;  a = ap.a ; z = ap.z ; invΔ = ap.invΔ ;  ρ = ap.ρ ;     z = ap.z ; a = ap.a ;
    C = aa.C  ; A = aa.A ; u = aa.u ; B = aa.B ;
    an = length(a)
    zn = length(z)
    invda = (an-1)/(a[end]-a[1])
    invγ = 1/γ
    inv1γ = 1.0/(1.0-γ)
    Cvals = nonzeros(C)
    Avals = nonzeros(A)
    copy!(Avals, Cvals)
    ij = zero(Int)
    @inbounds for zi in 1:zn, ai in 1:an
        ij += 1
        krange = nzrange(A, ij)
        rows = rowvals(A)[krange]
        if ai < an
            ∂V = (as.V[ij+1] - as.V[ij]) * invda
        else
            # state constraint boundary condition
            ∂V = (as.w * z[zi] + as.r * a[end])^(-γ)
        end
        saving = as.w * z[zi] + as.r * a[ai] - ∂V^(-invγ)
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
                ∂V = (as.V[ij] - as.V[ij-1]) * invda
            else
                # state constraint boundary condition
                ∂V =  (as.w * z[zi] + as.r * a[1])^(-γ) 
            end
            saving = as.w * z[zi] + as.r * a[ai] - ∂V^(-invγ)
            if saving < 0 # case of negative drift
                if ai > 1
                    current = saving * invda 
                    index = searchsortedfirst(rows, ij - 1)
                    Avals[krange[index]] -= current
                    index = searchsortedfirst(rows, ij)
                    Avals[krange[index]] += current
                end
            else
                ∂V = (as.w * z[zi] + as.r * a[ai])^(-γ)
            end
        end

        # update u
        u[ij] = inv1γ * ∂V^(1-invγ) + invΔ * as.V[ij]
    end
end

function update_B!(ap::AiyagariProblem, aa::AiyagariArrays, as::AiyagariSolution)
    γ = ap.γ  ;  a = ap.a ; z = ap.z ; invΔ = ap.invΔ ;  ρ = ap.ρ ;     z = ap.z ; a = ap.a ;
    C = aa.C  ; A = aa.A ; u = aa.u ; B = aa.B ;
    an = length(a)
    zn = length(z)
    Bvals = nonzeros(B) 
    Brows = rowvals(B)
    Avals = nonzeros(A)
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
end

function solve_hjb!(ap::AiyagariProblem, 
                    aa::AiyagariArrays,
                    as::AiyagariSolution; 
                    maxit::Int = 100, 
                    crit::Float64 = 1e-6, 
                    verbose::Bool = true)
    # update price vectors

    for iter in 1:maxit
        # update newV
        update_Au!(ap, aa, as)
        update_B!(ap, aa, as)
        newV = aa.B' \ aa.u
        # check convergence
        distance = chebyshev(newV, as.V)
        if distance < crit
            if verbose
                println("hjb solved : $(iter) iterations")
            end
            break
        else
            # update V using newV
            (newV, as.V) = (as.V, newV)
        end
    end
end


##############################################################################
##
## Solve fokker-planck equation
##
##############################################################################

# solve fokker-planck
function solve_fp!(ap::AiyagariProblem, aa::AiyagariArrays, as::AiyagariSolution)
    A = aa.A 
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
    as.g = A \ aa.b
    scale!(as.g, 1/sum(as.g))
end


##############################################################################
##
## Solve equilibrium by iterating on K
##
##############################################################################
function update_rw!(ap::AiyagariProblem, as::AiyagariSolution)
    as.r = exp(ap.π) * ap.α * as.K^(ap.α-1) - ap.δ    
    as.w = exp(ap.π) * (1-ap.α) * as.K^ap.α   
end

# compute agregate capital
function sum_capital{T}(a::Vector{Float64}, g::Vector{T})
    numerator = zero(Float64)
    denominator = zero(Float64)
    an = length(a)
    zn = div(length(g), length(a))
    ij = zero(Int)
    @inbounds for zi in 1:zn
        @simd for ai in 1:an
            ij += 1
            current = g[ij]
            numerator += current * a[ai] 
            denominator += current
        end
    end
    return numerator / denominator
end

# solve
function solve!(ap::AiyagariProblem, aa::AiyagariArrays, as::AiyagariSolution; 
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
        update_rw!(ap, as)

        # solve hamilton-jacobi-bellman equation 
        solve_hjb!(ap, aa, as, maxit = maxit, crit = crit, verbose = verbose)

        # solve fokker-planck equation
        solve_fp!(ap, aa, as)

        # Update agregate capital
        newK = sum_capital(ap.a, as.g)
        if verbose
            @show newK
        end
        if abs(as.K - newK) < critK
            break
        else
            # relaxation algorithm (to ensure convergence)
            as.K = relax * as.K + (1 - relax) * newK  
        end
    end
    return as
end

function solve(ap::AiyagariProblem;
    maxit::Int  = 100,      # maximum number of iterations in the HznB loop
    maxitK::Int = 100,      # maximum number of iterations in the K loop
    crit::Float64 = 1e-6,   # criterion HznB loop
    critK::Float64 = 1e-5,  # criterion K loop
    relax::Float64 = 0.99,  # relaxation parameter 
    verbose = true)
    as = AiyagariSolution(ap)
    aa = AiyagariArrays(ap, as)
    solve!(ap, aa, as; maxit = maxit, maxitK = maxitK, crit = crit, critK = critK, relax = relax, verbose = verbose)
end

# as = AiyagariSolution(ap)
# HJBFiniteDifference.solve_hjb!(ap, as)
