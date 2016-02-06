module HJBFiniteDifference

import Distributions: Gamma, Normal
import StatsBase: denserank
import DataFrames: DataFrame
import Gadfly: plot, Geom
import NLsolve: nlsolve
import Distances: chebyshev
using ODE
using CompEcon

##############################################################################
##
## Load files
##
##############################################################################
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
export solve!,
solve,
AiyagariProblem,
AiyagariMethod,
AiyagariFD,
AiyagariSimple,
DynamicAiyagariProblem,
BansalYaronProblem,
plot_ll

end