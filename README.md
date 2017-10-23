

- The package solves the Aiyagari model following Achdou, Han, Lasry, Lions and Moll (2015) "Heterogeneous Agent Models in Continuous Time"
	```julia
	using AiyagariContinuousTime
	ap = AiyagariProblem(π = 0.0);
	as = solve(ap)
	```

- To install, 
	```julia
	Pkg.clone("https://github.com/QuantEcon/Gensys.jl")
	Pkg.clone("https://github.com/matthieugomez/AiyagariContinuousTime.jl")
	```

# Bibliography
[Heterogeneous Agent Models in Continuous Time](http://www.princeton.edu/~moll/HACTproject.htm) by Ben Moll.