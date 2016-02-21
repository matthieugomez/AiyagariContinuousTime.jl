##############################################################################
##
## Gensys solver adapted from phactsolver.m
##
##############################################################################

function gensys(Γ0, Γ1, c, Ψ, Π; clean = true, continuous = true, check_existence = true, check_uniqueness = true)

    if clean
        @sprintf "Converting to Reduced Form"
        redundant = (maxabs(Γ0, 2) .== 0) & (maxabs(Ψ, 2) .== 0)
        base = nullspace(Γ1[redundant, :])
        Γ0 = lufact!(At_mul_B(base, Γ0 * base))
        try
            Γ1 = Γ0 \ At_mul_B(base, Γ1 * base)
            Ψ = Γ0 \ At_mul_B(base, Ψ)
            Π = Γ0 \ At_mul_B(base, Π)
            c = Γ0 \ At_mul_B(base, c)
        catch
            error("Wrong Form. Try running Gensys")
        end
    else
        Γ1 = Γ0 \ Γ1
    end
    n = size(Γ1, 1)

    # Schur Decomposition
    Γ1 = schurfact!(Γ1)
    if continuous
        select = real(Γ1[:values]) .< 0
    else
        select = abs(Γ1[:values]) .< 1
    end
    ordschur!(Γ1, select)
    n1 = sum(select)
    Γ1vectors = Γ1[:vectors]
    G1 = real(
        A_mul_Bt(Γ1vectors * Γ1[:Schur] * diagm(vcat(ones(n1), zeros(n - n1))), Γ1vectors))

    # thin svd
    u2 = Γ1vectors[:, (n1 + 1):n]
    etawt = svdfact!(At_mul_B(u2, Π))
    ueta, deta, veta  = etawt[:U], etawt[:S], etawt[:V]
    impact = real(-Π * veta * (diagm(deta) \ ueta') * At_mul_B(u2, Ψ) + Ψ)

    # check existence
    if check_existence
        temp = svdfact!(At_mul_B(u2, Ψ))
        uz, dz, vz =  temp[:U], temp[:S], temp[:V]
        existence = vecnorm(uz - ueta * At_mul_B(ueta, uz), 2) < (sqrt(eps()) * 10 * n)
    end

    # check uniqueness
    if check_uniqueness
        u1 = Γ1vectors[:, 1:n1]
        temp = svdfact!(At_mul_B(u1, Π))
        dont, deta1, veta1 = temp[:U], temp[:S], temp[:V]
        uniqueness = vecnorm(veta1 - veta * At_mul_B(veta, veta1), 2) < (sqrt(eps()) * 10  *  n)
    end

    if clean
        G1 = base * A_mul_Bt(G1, base)
        impact = base * impact
    end

    return G1, impact, existence, uniqueness
end