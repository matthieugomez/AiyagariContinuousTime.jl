

- The package solves the Aiyagari model following Achdou, Han, Lasry, Lions and Moll (2015) "Heterogeneous Agent Models in Continuous Time"
	```julia
	using AiyagariContinuousTime
	# solve a static equilibrium
	ap = AiyagariProblem(π = 0.0);
	as = solve(ap)
	```
- The package solves a *untested* dynamic version of the Aiyagri model following Ahn, Kaplan, Moll, Winberry (2016) "How to solve heterogeneous agent models in continuous time, with aggregate shocks"

	```julia
	using AiyagariContinuousTime, Gadfly
	ap = AiyagariProblem(π = 0.0);
	ρπ = 0.95
	σπ = 0.007
	as = solve(ap, ρπ, σπ)
	time, V, g, K, r, w = simulate(ap, as)
	plot(x = time, y = r, Geom.line,  Guide.xlabel("Years"), Guide.ylabel("Percentage points"), Guide.title("Interest Rate to Aggregate Productivity Shock"))
	```
	![aiyagari](https://cdn.rawgit.com/matthieugomez/AiyagariContinuousTime.jl/master/img/aiyagaridynamic.svg)

- To install, 
	```julia
	Pkg.clone("https://github.com/QuantEcon/Gensys.jl")
	Pkg.clone("https://github.com/matthieugomez/AiyagariContinuousTime.jl")
	```

# Bibliography
[Heterogeneous Agent Models in Continuous Time](http://www.princeton.edu/~moll/HACTproject.htm) by Ben Moll.