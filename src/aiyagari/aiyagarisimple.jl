##############################################################################
##
## Type
##
##############################################################################

type AiyagariSimple <:AiyagariMethod
    ap::AiyagariProblem 

    # discretization
    z::Vector{Float64}       # productivity vector
    invdz::Float64
    a::Vector{Float64}       # wealth vector
    invda::Float64

    # solution
    K::Float64               # aggregate capital 
    r::Float64
    w::Float64
    V::Vector{Float64}       # value function
    gg::Vector{Float64}      # distribution
end


function AiyagariSimple(ap::AiyagariProblem;  
                     zn::Int = 30, 
                     zmin::Float64 = 0.5, 
                     zmax::Float64 = 1.5, 
                     an::Int = 100, 
                     amin::Float64 = -1.0, 
                     amax::Float64 = 30.0, 
                     invΔ::Float64 = 1e-3,
                     K::Float64 = 3.8
                     )
    π = ap.π ; α = ap.α ; δ = ap.δ ; γ = ap.γ ; ρ = ap.ρ ;  σ2 = ap.σ2 ; θ = ap.θ

    a = collect(linspace(amin, amax, an))  
    invda = 1/ (a[2]-a[1])
    z = collect(linspace(zmin, zmax, zn))  
    invdz = 1/ (z[2]-z[1])

    # initial value function
    r = π * α * K^(α-1) - δ 
    w = π * (1-α) * K^α   
    V = Array(Float64, an*zn)
    ij = zero(Int)
    for zi in 1:zn, ai in 1:an
        ij += 1
        V[ij] = (w * z[zi] + r * a[ai])^(1-γ)/((1-γ)*ρ)
    end
    # storage arrays
    gg = ones(an*zn)
    scale!(gg, invda * invdz / sum(gg))
    AiyagariSimple(ap, z, invdz, a, invda, K, r, w, V, gg)
end



##############################################################################
##
## Solve hamilton-jacobi-bellman
##
##############################################################################
function V∂a(as, V, ai, zi)
    z = as.z ; a = as.a ; w = as.w ; r = as.r; γ = as.ap.γ; invda = as.invda ; invdz = as.invdz
    an = length(a)
    zn = length(z)
    ij = an * (zi - 1) + ai
    V∂ai = zero(eltype(V))
    if ai == 1
        # state constraint
        V∂ai = (z[zi] * w + a[ai] * r)^(-γ)
    elseif ai == an
        V∂ai = 0.5 * (3 * V[ij] - 4 * V[ij - 1] + V[ij - 2]) * invda
    else
        V∂ai = 0.5 * (V[ij + 1] - V[ij - 1]) * invda
    end
    return V∂ai
end

function V∂z(as, V, ai, zi)
    z = as.z ; a = as.a ; w = as.w ; r = as.r; γ = as.ap.γ; invda = as.invda ; invdz = as.invdz
    an = length(a)
    zn = length(z)
    ij = an * (zi - 1) + ai
    # V∂z = 0 at the border
    V∂zi = zero(eltype(V))
    if (zi > 1) & (zi < zn)
        V∂zi = 0.5 * (V[ij + an] - V[ij - an]) * invdz
    end
    return V∂zi
end

function V∂2z(as, V, ai, zi)
    z = as.z ; a = as.a ; w = as.w ; r = as.r; γ = as.ap.γ; invda = as.invda ; invdz = as.invdz
    an = length(a)
    zn = length(z)
    ij = an * (zi - 1) + ai
    V∂2zi = zero(eltype(V))
    if zi == 1
       V∂2zi = (2 * V[ij + an] - 2 * V[ij]) * invdz^2
    elseif zi == zn
       V∂2zi = (2 * V[ij - an] - 2 * V[ij]) * invdz^2
    else
       V∂2zi = (V[ij + an] + V[ij - an] - 2 * V[ij]) * invdz^2
    end
    return V∂2zi
end


function F!(as::AiyagariSimple, V, Vdot)
    π = as.ap.π ; α = as.ap.α ; δ = as.ap.δ ; γ = as.ap.γ ; ρ = as.ap.ρ ;  σ2 = as.ap.σ2 ; θ = as.ap.θ
    z = as.z ; a = as.a ; w = as.w ; r = as.r
    an = length(a)
    zn = length(z)
    invγ = 1/γ
    inv1γ = 1.0/(1.0-γ)
    ij = zero(Int)
    @inbounds for zi in 1:zn, ai in 1:an
        ij += 1     
        Vi = V[ij]
        V∂ai = V∂a(as, V, ai, zi)
        V∂zi = V∂z(as, V, ai, zi)
        V∂2zi =  V∂2z(as, V, ai, zi)
        Vdot[ij] = inv1γ * V∂ai^(1-invγ) + V∂ai * (w * z[zi] + r * a[ai] - V∂ai^(-invγ)) + θ * (1.0 - z[zi]) * V∂zi + 0.5 * σ2 * V∂2zi - ρ * Vi
    end
end

function solve_hjb!(as::AiyagariSimple;  
                    maxit::Int = 100, 
                    crit::Float64 = 1e-6, 
                    verbose::Bool = true)
    # update price vectors
    as.r = as.ap.π * as.ap.α * as.K^(as.ap.α-1) - as.ap.δ    
    as.w = as.ap.π * (1-as.ap.α) * as.K^as.ap.α     
    out = nlsolve((V, Vdot) -> F!(as, V, Vdot), as.V, autodiff = true, xtol = 1e-10, ftol = 1e-10)
    as.V = out.zero
end


##############################################################################
##
## Solve fokker-planck equation
##
##############################################################################

# solve fokker-planck
function G!(as::AiyagariSimple, y, ydot)
    π = as.ap.π ; α = as.ap.α ; δ = as.ap.δ ; γ = as.ap.γ ; ρ = as.ap.ρ ;  σ2 = as.ap.σ2 ; θ = as.ap.θ
    z = as.z ; a = as.a ; w = as.w ; r = as.r ; invda = as.invda ; invdz = as.invdz
    an = length(a)
    zn = length(z)
    invγ = 1/γ
    ij = zero(Int)
    for zi in 1:zn, ai in 1:an
        ij += 1     
        if ij == 2
            ydot[ij] = y[ij] - 0.1
        else
            g∂ai = zero(eltype(y))
            if ai < an
                g∂ai += (w * z[zi] + r * a[ai+1] - V∂a(as, as.V, ai + 1, zi)^(-invγ)) * 0.5 * y[ij + 1] * invda
            end
            if ai > 1
                g∂ai += - (w * z[zi] + r * a[ai-1] - V∂a(as, as.V, ai - 1, zi)^(-invγ)) * 0.5 * y[ij - 1] * invda
            end
            g∂zi = zero(eltype(y))
            if zi < zn
                g∂zi += θ * (1.0 - z[zi + 1]) * 0.5 * y[ij + an] * invdz
            end
            if zi > 1
                g∂zi += - θ * (1.0 - z[zi - 1]) * 0.5 * y[ij - an] * invdz
            end
            g∂2zi = - 0.5 * σ2 * 2 * y[ij] * invdz^2
            if zi < zn
                g∂2zi += 0.5 * σ2 * y[ij + an] * invdz^2
            end
            if zi > 1
                g∂2zi += 0.5 * σ2 * y[ij - an] * invdz^2
            end
            ydot[ij] = - g∂ai - g∂zi + g∂2zi
        end
    end
    return ydot
end

function solve_fp!(as::AiyagariSimple)
    out = nlsolve((g, gdot) -> G!(as, g, gdot), deepcopy(as.gg), autodiff = true)
    @show out
    as.gg = out.zero
end
    







