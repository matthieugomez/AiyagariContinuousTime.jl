#=
Solve the Huggett economy from

Achdou, Y., Han, J., Lasry, J.-M., Lions, P.-L., & Moll, B. (2015).
Heterogeneous Agent Models in Continuous Time.
Retrieved from http://www.princeton.edu/
=#

# ---------------- #
# Model Parameters #
# ---------------- #

immutable HuggettProblem
    # discount and interest rates
    ρ::Float64
    r::Float64

    # CRRA parameter
    γ::Float64

    # labor income parameters
    λ1::Float64
    λ2::Float64
    z1::Float64
    z2::Float64

    # grid on assets. lower bound is borrowing constraint
    agrid::Vector{Float64}
    Δa::Float64
end

function HuggettProblem(;ρ=0.05, r=0.03, γ=2.0, λ1=0.02, λ2=0.03, z1=0.1, z2=0.2,
                         abar=-0.02, amax=2.0, Na=500)
    agrid = collect(linspace(abar, amax, Na))
    Δa = agrid[2] - agrid[1]
    HuggettProblem(ρ, r, γ, λ1, λ2, z1, z2, agrid, Δa)
end

function initial_vf!(m::HuggettProblem, v1::AbstractVector, v2::AbstractVector)
    for i in 1:length(m.agrid)
        v1[i] = (m.z1 + m.r*m.agrid[i]).^(1-m.γ)/(m.ρ*(1-m.γ))
        v2[i] = (m.z2 + m.r*m.agrid[i]).^(1-m.γ)/(m.ρ*(1-m.γ))
    end
    v1, v2
end

function explicit_solve(m::HuggettProblem, maxit::Int=20_000,
                        tol::Float64=1e-6)
    # extract parameters
    Na = length(m.agrid)
    m1_γ = -1/m.γ
    z = [m.z1, m.z2]
    λ = [m.λ1, m.λ2]

    # to be used later
    css = z' .+ m.r*m.agrid
    dvss = css.^(-m.γ)
    dvfN = zeros(z)
    dvb1 = (z + m.r*m.agrid[1]).^(-m.γ)

    # initial guess at vf
    v = Array(Float64, length(m.agrid), 2)
    initial_vf!(m, sub(v, :, 1), sub(v, :, 2))
    old_v = copy(v)

    start_t = time()

    for it in 1:maxit
        # reset err
        err = 0.0
        copy!(old_v, v)

        @inbounds for j in 1:2
            not_j = j == 1 ? 2 : 1
            for i in 1:Na
                dvf = i == Na ? dvfN[j] : (old_v[i+1, j] - old_v[i, j]) / m.Δa
                dvb = i == 1  ? dvb1[j] : (old_v[i, j] - old_v[i-1, j]) / m.Δa

                # @show it, j, i, dvf, dvb

                # use eqn 9 to update c and use budget constraint for s
                cb = dvb^(m1_γ)
                sb = z[j] + m.r*m.agrid[i] - cb

                cf = dvf^(m1_γ)
                sf = z[j] + m.r*m.agrid[i] - cf

                (dvij, cij, sij) = sf > 0 ? (dvf, cf, sf) :
                                   sb < 0 ? (dvb, cb, sb) :
                                            (dvss[i, j], css[i, j], 0.0)

                # compute rhs of 10 and update guess for v
                rhs = cij^(1-m.γ)/(1-m.γ) +
                      dvij*sij +
                      λ[j]*(old_v[i, not_j] - old_v[i, j]) -
                      m.ρ*old_v[i, j]
                v[i, j] = Δ*rhs + v[i, j]

                # track err as we go
                err = max(err, abs(rhs))
            end
        end

        # if it % 100 == 0
        #     @printf "%-8i%-12.5e%12.5f\n" it err (time() - start_t)
        # end

        if err < tol
            println("converged in $it steps in $(time() -start_t) seconds")
            return v
        end
    end

    error("maxit exceeded")
end

function implicit_solve(m::HuggettProblem, Δ::Float64=1000.0, maxit=10_000,
                        tol=1e-6)
    # extract parameters
    Na = length(m.agrid)
    m1_γ = -1/m.γ
    z = [m.z1, m.z2]
    λ = [m.λ1, m.λ2]

    # to be used later
    income = z' .+ m.r*m.agrid
    css = income
    dvss = css.^(-m.γ)
    dvfN = zeros(z)
    dvb1 = (z + m.r*m.agrid[1]).^(-m.γ)

    # initial guess at vf
    v = Array(Float64, length(m.agrid), 2)
    initial_vf!(m, sub(v, :, 1), sub(v, :, 2))
    old_v = copy(v)

    # lower, mid, and upper diagonal for Tridiagonal part of A
    dl = zeros(Float64, 2Na-1)
    d = zeros(Float64, 2Na)
    du = zeros(Float64, 2Na-1)
    b = similar(d)

    A_λ = [-speye(Na)*λ[1] speye(Na)*λ[1];speye(Na)*λ[2] -speye(Na)*λ[2]]

    start_t = time()

    @inbounds for it in 1:maxit
        copy!(old_v, v)
        for j in 1:length(z)
            for i in 1:Na
                ind = sub2ind((Na, 2), i, j)
                dvf = i == Na ? dvfN[j] : (old_v[i+1, j] - old_v[i, j]) / m.Δa
                dvb = i == 1  ? dvb1[j] : (old_v[i, j] - old_v[i-1, j]) / m.Δa

                # use eqn 9 to update c and use budget constraint for s
                cb = dvb^(m1_γ)
                sb = income[i, j] - cb

                cf = dvf^(m1_γ)
                sf = income[i, j] - cf

                uij = (sf > 0 ? cf :
                       sb < 0 ? cb :
                                css[i, j])^(1-m.γ)/(1-m.γ)

                # fill diagonal with y{i,j}
                d[ind] = (- max(sf, 0.0) + min(sb, 0.0))/m.Δa

                # fill in lower and upper diagonals.
                # when i is 1 we don't use a lower diagonal.
                # When i == Na we don't use an upper
                # otherwise use both
                if i == 1
                    du[ind] = max(sf, 0.0)/m.Δa
                elseif i == Na
                    dl[ind-1] = - min(sb, 0.0)/m.Δa
                else
                    du[ind] = max(sf, 0.0)/m.Δa
                    dl[ind-1] = - min(sb, 0.0)/m.Δa
                end

                b[ind] = uij + old_v[ind]/Δ

            end

        end

        A = sparse(Tridiagonal(dl, d, du)) + A_λ
        B = (1/Δ + m.ρ)I - A

        copy!(v, B\b)
        err = chebyshev(v, old_v)

        if err < tol
            return v
        end

        @printf "%-8i%-12.5e%12.5f\n" it err (time() - start_t)

    end

    error("maxit exceeded")
end


function solve(m::HuggettProblem; method::Symbol=:implicit, maxit::Int=20_000,
               tol::Float64=1e-6, Δ::Float64=1000.0)

    if method == :implicit
        implicit_solve(m, Δ, maxit, tol)
    else
        explicit_solve(m, maxit, tol)
    end
end
