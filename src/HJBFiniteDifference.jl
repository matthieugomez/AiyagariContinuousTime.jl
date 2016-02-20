module HJBFiniteDifference

import Distributions: Gamma, Normal
import StatsBase: denserank
import DataFrames: DataFrame
import Gadfly: plot, Geom
import NLsolve: nlsolve
import Distances: chebyshev
using ForwardDiff
using ODE
using CompEcon

##############################################################################
##
## Load files
##
##############################################################################
include("aiyagari/aiyagari.jl")
include("aiyagari/dynamicaiyagari.jl")
include("aiyagari/phact_solver.jl")

include("bansalyaron/bansalyaronproblem.jl")
include("bansalyaron/finitedifferences.jl")
include("bansalyaron/spectralmethod.jl")



##############################################################################
##
## Exported methods and types 
##
##############################################################################
export solve,
AiyagariProblem,
AiyagariArrays,
AiyagariSolution,
AiyagariMethod,
BansalYaronProblem,
plot_ll
#phact_solver

end