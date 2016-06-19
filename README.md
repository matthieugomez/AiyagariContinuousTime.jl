# Install
```julia
Pkg.clone("https://github.com/matthieugomez/HJBFiniteDifference.jl")
```

# Kolmogorov Forward
The package solves the kolmogorov forward equation on a grid, i.e. the stationary distribution `g` that solves

`0 = -∂(μ g) + 0.5 * ∂^2(σ^2 g) + δ (ψ - 1)`

The function outputs a vector `g Δx`  which sums to one. The grid `x` is potentially non-uniform. If `σ` is not null at the boundary of the grids, the process is assumed to be reflected.

For instance, let us plot the stationary distribution of a brownian motion with negative drift (a pareto distribution):
```julia
using HJBFiniteDifference
x = logspace(-2, 5, 100)
μ = -0.001 .* x
σ = 0.2 .* x
g = kolmogorovforward(x, μ, σ)
cumulativeg = cumsum(g)
using Gadfly
plot(x = log(x), y = log(1 .- cumulativeg), Geom.line, Guide.xlabel("log-x"), Guide.ylabel("log 1-cdf"))
```
![kolmogorov](https://cdn.rawgit.com/matthieugomez/HJBFiniteDifference.jl/master/img/kolmogorov.svg)



# Aiyagari
- The package solves the Aiyagari model following Achdou, Han, Lasry, Lions and Moll (2015) "Heterogeneous Agent Models in Continuous Time"
```julia
using HJBFiniteDifference
# solve a static equilibrium
ap = AiyagariProblem(π = 0.0);
as = solve(ap)
```
- The package solves a dynamic version of the Aiyagri model following [Ahn, Kaplan, Moll, Winberry](How to solve heterogeneous agent models in continuous time, with aggregate shocks)

```julia
using HJBFiniteDifference, Gadfly
ap = AiyagariProblem(π = 0.0);
ρπ = 0.95
σπ = 0.007
as = solve(ap, ρπ, σπ)
time, V, g, K, r, w = simulate(ap, as)
plot(x = time, y = r, Geom.line,  Guide.xlabel("Years"), Guide.ylabel("Percentage points"), Guide.title("Interest Rate to Aggregate Productivity Shock"))
```
![aiyagari](https://cdn.rawgit.com/matthieugomez/HJBFiniteDifference.jl/master/img/aiyagaridynamic.svg)


# Bansal Yaron

The package solves the PDE associated with the long run risk model of Bansal-Yaron (2004). This long run risk model is generally solved by log-linearization (i.e. assuming that the price dividend is log linear in state variables). Solving directly the PDE shows that the price-dividend actually displays substantial non linearity wrt volatility. The choice of parameters follows Bansal-Kiku-Yaron (2007). Explanations for the solution method are available [here](https://github.com/matthieugomez/HJBFiniteDifference.jl/blob/master/src/bansalyaron/bansalyaron.pdf).


```julia
using HJBFiniteDifference, Gadfly
byp = BansalYaronProblem(ρ = -log(0.9989), γ = 7.5, ψ = 1.5, νD = 0.0072, νμ = 0.038 * 0.0072, νσ = 0.00000283 / 0.0072^2, κμ = -log(0.975), κσ = -log(0.999))

# Finite Differences (Nonlinear solver)
solution = solve(byp, method = :nl)
plot(byp, solution, :s2)
plot(byp, solution, :m)

# Finite Differences (ODE solver)
solution = solve(byp, method = :ode)
plot(byp, solution, :s2)
plot(byp, solution, :m)

# Spectral Method (Chebyshev polynomials)
solution = solve(byp, method = :spectral)
plot(byp, solution, :s2)
plot(byp, solution, :m)
```
![bansalyaron](https://cdn.rawgit.com/matthieugomez/HJBFiniteDifference.jl/master/img/byp_m.svg)
![bansalyaron](https://cdn.rawgit.com/matthieugomez/HJBFiniteDifference.jl/master/img/byp_vol.svg)


# Bibliography
Three excellent resources to learn about finite difference schemes and their applications to HJB equations:
- [Numerical analysis of partial differential equations arising in finance and stochastic control](http://www.cmap.polytechnique.fr/%7Ebonnans/notes/edpfin/edpfin.html) by Frédéric Bonnans.
- [An introduction to Finite Difference methods for
PDEs in Finance](https://www.fields.utoronto.ca/programs/scientific/09-10/finance/courses/tourin.pdf)  by Agnès Tourin 
-  [Heterogeneous Agent Models in Continuous Time](http://www.princeton.edu/~moll/HACTproject.htm) by Ben Moll.