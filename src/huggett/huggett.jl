#=
Solve the Huggett economy from

Achdou, Y., Han, J., Lasry, J.-M., Lions, P.-L., & Moll, B. (2015).
Heterogeneous Agent Models in Continuous Time.
Retrieved from http://www.princeton.edu/
=#

# ---------------- #
# Model Parameters #
# ---------------- #

type HuggettProblem
    # discount and interest rates
    ρ::Float64
    r::Float64

    # CRRA parameter
    γ::Float64

    # labor income parameters
    λ::Vector{Float64}
    z::Vector{Float64}

    # grid on assets. lower bound is borrowing constraint
    agrid::Vector{Float64}
    Δa::Float64
end

function HuggettProblem(;ρ=0.05, r=0.03, γ=2.0, λ=[1.2, 1.2], z=[0.1, 0.2],
                         abar=-0.15, amax=5.0, Na=1000)
    agrid = collect(linspace(abar, amax, Na))
    Δa = agrid[2] - agrid[1]
    HuggettProblem(ρ, r, γ, λ, z, agrid, Δa)
end

Base.copy(m::HuggettProblem) =
    HuggettProblem([copy(getfield(m, k)) for k in fieldnames(m)]...)

function initial_vf!(m::HuggettProblem, v)
    size(v) != (length(m.agrid), length(m.z)) && error("v wrong shape")
    for j in 1:length(m.z), i in 1:length(m.agrid)
        v[i, j] = (m.z[j] + m.r*m.agrid[i]).^(1-m.γ)/(m.ρ*(1-m.γ))
        v[i, j] = (m.z[j] + m.r*m.agrid[i]).^(1-m.γ)/(m.ρ*(1-m.γ))
    end
    v
end

initial_vf(m) = initial_vf!(m, Array(Float64, length(m.agrid), length(m.z)))

function explicit_solve_hjb(m::HuggettProblem, v=initial_vf(m),
                            maxit::Int=20_000, tol::Float64=1e-6, verbose=false)
    # extract parameters
    Na, z, λ = length(m.agrid), m.z, m.λ
    m1_γ = -1/m.γ

    # to be used later
    css = z' .+ m.r*m.agrid
    dvss = css.^(-m.γ)
    dvfN = zeros(z)
    dvb1 = (z + m.r*m.agrid[1]).^(-m.γ)

    # initial guess at vf
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
            if verbose
                println("converged in $it steps in $(time() -start_t) seconds")
            end
            # a hack to get the A matrix from the implicit method so we can
            # solve FP below.
            return v, implicit_solve_hjb(m, v)[2]
        end
    end

    error("maxit exceeded")
end

function implicit_solve_hjb(m::HuggettProblem, v=initial_vf(m),
                            Δ::Float64=1000.0, maxit=10_000, tol=1e-6,
                            verbose=false)
    # extract parameters
    Na, z, λ = length(m.agrid), m.z, m.λ
    m1_γ = -1/m.γ

    # to be used later
    income = z' .+ m.r*m.agrid
    css = income
    dvss = css.^(-m.γ)
    dvfN = zeros(z)
    dvb1 = (z + m.r*m.agrid[1]).^(-m.γ)

    # initial guess at vf
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
            return v, A
        end

        if verbose
            @printf "%-8i%-12.5e%12.5f\n" it err (time() - start_t)
        end

    end

    error("maxit exceeded")
end

function solve_kolmogorov_forward(m::HuggettProblem, A::AbstractSparseMatrix)
    Na = length(m.agrid)
    AT = A'
    rhs = zeros(2Na)

    # fix one eigenvalue, otherwise matrix is singular
    fixer = 1
    rhs[fixer] = 0.1
    AT[fixer, :] = 0.0
    AT[fixer, fixer] = 1.0

    # normalize
    gg = AT\rhs
    gg ./= (sum(gg)*m.Δa)

    # put into columns
    g = [gg[1:Na] gg[Na+1:2Na]]

    # check that we actually have what we think we do
    @assert sum(sum(g, 1) .* m.Δa) - 1.0 < 1e-10

    g
end


function solve_hjb(m::HuggettProblem; method::Symbol=:implicit, maxit::Int=20_000,
                   tol::Float64=1e-6, Δ::Float64=1000.0, verbose=false)

    if method == :implicit
        implicit_solve_hjb(m, initial_vf(m), Δ, maxit, tol, verbose)
    else
        explicit_solve_hjb(m, initial_vf(m), maxit, tol, verbose)
    end
end

function solve_partial_eq(m::HuggettProblem; method::Symbol=:implicit,
                          maxit::Int=20_000, tol::Float64=1e-6,
                          Δ::Float64=1000.0, verbose=false)
    v, A = solve_hjb(m; method=method, maxit=maxit, tol=tol, Δ=Δ,
                     verbose=verbose)
    g = solve_kolmogorov_forward(m, A)
    v, g, A
end

aggregate_savings(m::HuggettProblem, v, g) = sum(sum(g'm.agrid)) * m.Δa

function solve_general_eq(m::HuggettProblem, rl=0.01, ru=0.04;
                          verbose::Bool=false)
    rstar = brent(rl, ru) do r
        m_r = copy(m); m_r.r = r
        v, g, A = solve_partial_eq(m_r; verbose=verbose)
        aggregate_savings(m, v, g)
    end

    # build model with this interest rate
    m_out = copy(m); m_out.r = rstar

    # build partial equilibrium solution
    v, g, A = solve_partial_eq(m_out; verbose=verbose)

    # reutrn model and solution
    m, v, g, A

end
