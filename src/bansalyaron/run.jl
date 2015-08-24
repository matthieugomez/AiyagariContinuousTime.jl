import DataFrames: DataFrame
import Gadfly: plot, Geom
import NLsolve: nlsolve, DifferentiableSparseMultivariateFunction
abstract BansalYaronProblem
include("/Users/Matthieu/Dropbox/Personal/Research/Code/Julia/finitedifference/bansalyaron/utils.jl")
include("/Users/Matthieu/Dropbox/Personal/Research/Code/Julia/finitedifference/bansalyaron/bansalyaronL.jl")
include("/Users/Matthieu/Dropbox/Personal/Research/Code/Julia/finitedifference/bansalyaron/bansalyaronNL.jl")





#
# Bansal Yaron (2004)
#

# linear scheme (runs in 0.06s)
byp =	BansalYaronProblemL(γ = 7.5, ψ = 1.5)
@time solve_hmj!(byp)
plot(byp, :s2)
plot(byp, :m)

# non linear scheme (runs in 0.6s)
byp =	BansalYaronProblemNL(γ = 7.5, ψ = 1.5)
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


