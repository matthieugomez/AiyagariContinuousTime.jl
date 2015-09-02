# Install
```julia
Pkg.clone("https://github.com/matthieugomez/HJBviaPDE.jl")
```

# Aiyagari
Following Achdou, Han, Lasry, Lions and Moll (2015) "Heterogeneous Agent Models in Continuous Time"

```julia
# solve a static equilibrium
using HJBviaPDE, Gadfly
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

![aiyagari](https://cdn.rawgit.com/matthieugomez/HJBviaPDE.jl/master/img/aiyagari.svg)


# Bansal Yaron

Parameters from Bansal Yaron Kiku (2007)

```julia
using HJBviaPDE, Gadfly

# NL scheme solved using Newton method
byp = BansalYaronProblemNewton(ρ = -log(0.9989), γ = 7.5, ψ = 1.5, νD = 0.0072, νμ = 0.038 * 0.0072, νσ = 0.0000028 / 0.0072^2, κμ = -log(0.975), κσ = -log(0.999))
@time solve!(byp)
plot(byp, :s2)
plot(byp, :m)

# NL scheme solved using Powell Dog-leg method
byp = BansalYaronProblemPowell(ρ = -log(0.9989), γ = 7.5, ψ = 1.5, νD = 0.0072, νμ = 0.038 * 0.0072, νσ = 0.0000028 / 0.0072^2, κμ = -log(0.975), κσ = -log(0.999))
@time solve!(byp)
plot(byp, :s2)
plot(byp, :m)
```

![bansalyaron](https://cdn.rawgit.com/matthieugomez/HJBviaPDE.jl/master/img/bansalyaron.svg)


# Bibliography
Good resources on finite difference scheme to solve HJB numerically:
- *Numerical analysis of partial differential equations arising in finance and stochastic control* by Frédéric Bonnans. [Part1](http://www.cmap.polytechnique.fr/~bonnans/notes/edpfin/edpfin.html) and [Part2](http://www.cmap.polytechnique.fr/~bonnans/notes/edpfin/edpfinFVEL-dec2014.pdf) of the course.
- Ben Moll's page [Heterogeneous Agent Models in Continuous Time](http://www.princeton.edu/~moll/HACTproject.htm)