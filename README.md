# Install
```julia
Pkg.clone("https://github.com/QuantEcon/Gensys.jl")
Pkg.clone("https://github.com/matthieugomez/AiyagariContinuousTime.jl")
```

# Kolmogorov Forward
The package solves the kolmogorov forward equation on a grid, i.e. the stationary distribution `g` that solves

`0 = -∂(μ g) + 0.5 * ∂^2(σ^2 g) + δ (ψ - 1)`

The function outputs a vector `g Δx`  which sums to one. The grid `x` is potentially non-uniform. If `σ` is not null at the boundary of the grids, the process is assumed to be reflected.

```julia

# stationary distribution of Brownian motion
using HJBFiniteDifference
x = logspace(-4, 3, 100)
μ = zeros(x)
σ = 0.3 .* x
g = kolmogorovforward(x, μ, σ)
plot(x = log(x), y = log(1 .- cumsum(g)), Geom.line, Guide.xlabel("log-x"), Guide.ylabel("log 1-cdf"))

# stationary distribution of Ornstein-Uhlenbeck process
x = collect(linspace(0.0, 0.04, 10))
μ = [0.3 * (0.02 - x) for x in μ]
σ = [0.01 for x in μ]
g = kolmogorovforward(μ, μμ, σμ)
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
- The package solves a (untested) dynamic version of the Aiyagri model following Ahn, Kaplan, Moll, Winberry (2016) "How to solve heterogeneous agent models in continuous time, with aggregate shocks"

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


# Bibliography
[Heterogeneous Agent Models in Continuous Time](http://www.princeton.edu/~moll/HACTproject.htm) by Ben Moll.