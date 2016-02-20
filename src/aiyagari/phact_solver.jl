##############################################################################
##
## PHACT solver
##
##############################################################################

function phact_solver(g0, g1, c, psi, pi; clean = true, continuous = true, check_exist = true, check_uniq = true)
    realsmall = sqrt(eps()) * 10
    eu = [false, false]
    if clean
        @sprintf "Converting to Reduced Form"
        redundant = maximum(abs([g0 psi]), 2) .== 0
        base = nullspace(g1[redundant, :])
        g0 = At_mul_B(base, g0 * base)
        g1 = At_mul_B(base, g1 * base)
        g1 = g0 \ g1
        psi = g0 \ At_mul_B(base, psi)
        pi = g0 \ At_mul_B(base, pi)
        c = g0 \ At_mul_B(base, c)
    else
        g1 = g0 \ g1
    end
    n = size(g1, 1)

    # Schur Decomposition
    S = schurfact!(g1)
    if continuous
        select = real(S[:values]) .< 0
        ordschur!(S, select)
        g_eigs = S[:values]
        nunstab = sum(real(g_eigs) .> 0)
    else
        select = abs(S[:values]) .< 1
        ordschur!(S, select)
        g_eigs = abs(S[:values])
        nunstab = sum(real(g_eigs) .< 1)
    end

    u2 = S[:vectors][:, (n - nunstab + 1):n]'
    u1 = S[:vectors][:, 1:(n - nunstab)]'

    etawt = svdfact!(u2 * pi)
    ueta, deta, veta  = etawt[:U], etawt[:S], etawt[:V]
    md = length(deta)
    bigev = find(deta[1:md] .> realsmall)
    ueta = ueta[:, bigev]
    veta = veta[:, bigev]
    deta = deta[bigev]

    if check_exist
        zwt = svdfact!(u2 * psi)
        uz, dz, vz =  zwt[:U], zwt[:S], zwt[:V]
        md = length(dz)
        bigev = find(dz[1:md] .> realsmall)
        uz = uz[:, bigev]
        vz = vz[:, bigev]
        dz = dz[bigev, bigev]

        if isempty(bigev)
            eu[1] = true
        else
            eu[1] = vecnorm(uz - ueta * ueta' * uz, 2) < (realsmall * n)
        end

        if !eu[1]
            @sprintf "Solution does not exist"
        end
        impact = real(-pi * veta * (diagm(deta) \ ueta') * uz * dz * vz' + psi)
    else
        impact = real(-pi * veta * (diagm(deta) \ ueta') * u2 * psi + psi)
    end

    G1 = S[:vectors] * S[:Schur] * diagm(vcat(ones(n - nunstab), zeros(nunstab))) * S[:vectors]'
    G1 = real(G1)
        
    if check_uniq
        dont, deta1, veta1 = svd(u1 * pi)
        bigev = find(deta1 .> realsmall)
        veta1 = veta1[:, bigev]
        if isempty(veta1)
            eu[2] = true
        else
            eu[2] = vecnorm(veta1 - veta * veta' * veta1, 2) < (realsmall  *  n)
        end
    end

    if clean
        G1 = base * G1 * base'
        impact = base * impact
    end
    C = 1
    return G1, C, impact, eu
end