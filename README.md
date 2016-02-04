# Install
```julia
Pkg.clone("https://github.com/matthieugomez/HJBFiniteDifference.jl")
```

# Aiyagari
Following Achdou, Han, Lasry, Lions and Moll (2015) "Heterogeneous Agent Models in Continuous Time"

```julia
# solve a static equilibrium
using HJBFiniteDifference, Gadfly
π = 1.0
K = 3.8
ap = AiyagariProblem(π, K);
@time solve!(ap)

# solve a dynamic equilibrium
## construct an unexpected productivity shock π
dt = 0.5
π = Array(Float64, 400);
π[1] = 0.97
for t in 2:length(π)
    π[t] = 0.5 * 0.2 * (1.0 - π[t-1]) + π[t-1]
end
## solve along this shock
K = 3.8
apd = DynamicAiyagariProblem(π, K, dt = dt);
@time solve!(apd)
plot(x = collect(1:length(π)), y = apd.K, Geom.line)
```

![aiyagari](https://cdn.rawgit.com/matthieugomez/HJBFiniteDifference.jl/master/img/aiyagari.svg)


# Bansal Yaron

The package solves the PDE associated with the long run risk model of Bansal-Yaron (2004). This long run risk model is generally solved by log-linearization, which assumes that the price dividend is linear in state variables. Solving directly the PDE shows that the price-dividend actually displays substantial non linearity wrt volatility. The choice of parameters follows Bansal-Kiku-Yaron (2007).


```julia
using HJBFiniteDifference, Gadfly
byp = BansalYaronProblem(ρ = -log(0.9989), γ = 7.5, ψ = 1.5, νD = 0.0072, νμ = 0.038 * 0.0072, νσ = 0.0000028 / 0.0072^2, κμ = -log(0.975), κσ = -log(0.999))
solution = solve(byp, method = :spectral)

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

![bansalyaron](https://cdn.rawgit.com/matthieugomez/HJBFiniteDifference.jl/master/img/byp.svg)


# Bibliography
Two excellent resources to learn about finite difference schemes and their applications to HJB equations:
- [Numerical analysis of partial differential equations arising in finance and stochastic control](http://www.cmap.polytechnique.fr/%7Ebonnans/notes/edpfin/edpfin.html) by Frédéric Bonnans.
-  [Heterogeneous Agent Models in Continuous Time](http://www.princeton.edu/~moll/HACTproject.htm) by Ben Moll.