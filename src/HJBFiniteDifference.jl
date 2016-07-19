module HJBFiniteDifference

import Distributions: Gamma, Normal
import StatsBase: denserank
import DataFrames: DataFrame
import Gadfly: plot, Geom
import NLsolve: nlsolve
import Distances: chebyshev
using Gensys
using ForwardDiff
using ODE
using CompEcon
using Gensys

##############################################################################
##
## Load files
##
##############################################################################
include("kolmogorovforward/kolmogorov.jl")
include("aiyagari/aiyagari.jl")
include("aiyagari/dynamicaiyagari.jl")


include("bansalyaron/bansalyaronproblem.jl")
include("bansalyaron/finitedifferences.jl")
include("bansalyaron/spectralmethod.jl")



##############################################################################
##
## Exported methods and types 
##
##############################################################################
export solve,
kolmogorovforward,
simulate,
AiyagariProblem,
AiyagariArrays,
AiyagariSolution,
DynamicAiyagariSolution,
BansalYaronProblem,
plot_ll
#phact_solver

end