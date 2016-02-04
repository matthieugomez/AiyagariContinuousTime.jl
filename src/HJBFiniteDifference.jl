module HJBFiniteDifference


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
include("aiyagari/aiyagarisimple.jl")


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