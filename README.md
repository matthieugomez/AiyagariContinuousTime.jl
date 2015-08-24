This repository solves some economic models via PDE.


# Aiyagari
Following Achdou, Han, Lasry, Lions and Moll (2014)

```julia
# solve a static problem
π = 1.0
K = 3.8
ap = AiyagariProblem(π, K);
@time solve!(ap)


# solve a dynamic problem
N = 400 ;
dt = 0.5 ; 
ν = 1-0.8 ;
π = Array(Float64, N) ;
π[1]=.97
for t=1:(N-1)
    π[t+1] = dt * ν * (1 - π[t]) + π[t]
end
K = 3.8

# initializing computes K(t) stationary solution for π[t]
apd = DynamicAiyagariProblem(π, K, dt = dt);
solve!(apd)
#
## plot capital evolution
## install the package Gadfly with Pkg.add("Gadfly")
#using Gadfly
#plot(x = collect(1:N), y = apd.K, Geom.line)
```

# Bansal Yaron


```julia
#
# Bansal Yaron (2004)
#

# linear scheme (runs in 0.06s)
byp = BansalYaronProblemL(γ = 7.5, ψ = 1.5)
@time solve_hmj!(byp)
plot(byp, :s2)
plot(byp, :m)

# non linear scheme (runs in 0.6s)
byp = BansalYaronProblemNL(γ = 7.5, ψ = 1.5)
@time solve_hmj!(byp, method = :newton)
plot(byp, :s2)
plot(byp, :m)


#
# Bansal Yaron Kiku (2007)
#

# linear scheme (runs in 0.24s)
byp =BansalYaronProblemL(ρ = -log(0.9989), γ = 7.5, ψ = 1.5, νD = 0.0072, νμ = 0.038 * 0.0072, νσ = 0.0000028 / 0.0072^2, κμ = -log(0.975), κσ = -log(0.999))
@time solve_hmj!(byp)
plot(byp, :s2)
plot(byp, :m)

# non linear scheme (runs in 0.4s)
byp =BansalYaronProblemNL(ρ = -log(0.9989), γ = 7.5, ψ = 1.5, νD = 0.0072, νμ = 0.038 * 0.0072, νσ = 0.0000028 / 0.0072^2, κμ = -log(0.975), κσ = -log(0.999))
@time solve_hmj!(byp)
plot(byp, :s2)
plot(byp, :m)




# linear scheme (runs in 0.24s)
byp =BansalYaronProblemL(ρ = -log(0.9989), γ = 50, ψ = 1.5, νD = 0.0072, νμ = 0.038 * 0.0072, νσ = 0.0000028 / 0.0072^2, κμ = -log(0.975), κσ = -log(0.999))
@time solve_hmj!(byp)
plot(byp, :s2)
plot(byp, :m)

# non linear scheme (runs in 0.4s)
byp =BansalYaronProblemNL(ρ = -log(0.9989), γ = 50, ψ = 1.5, νD = 0.0072, νμ = 0.038 * 0.0072, νσ = 0.0000028 / 0.0072^2, κμ = -log(0.975), κσ = -log(0.999))
@time solve_hmj!(byp)
plot(byp, :s2)
plot(byp, :m)
```


