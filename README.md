
The package solves the Aiyagari model following Achdou, Han, Lasry, Lions and Moll (2015) "Heterogeneous Agent Models in Continuous Time"
	```julia
	using AiyagariContinuousTime
	ap = AiyagariProblem(π = 0.0);
	as = solve(ap)
	```

The package is no longer maintaned. Please go to the [EconPDEs](https://github.com/matthieugomez/EconPDEs.jl/tree/master/examples) repository for an updated version, that allows to solve a vast arrays of PDEs in Julia (including Consumption/Saving Models, Investment Models, and Asset Pricing Models).