module HJBFiniteDifference


import DataFrames: DataFrame
import Gadfly: plot, Geom
import NLsolve: nlsolve, DifferentiableSparseMultivariateFunction
import Distances: chebyshev


##############################################################################
##
## Exported methods and types 
##
##############################################################################
export solve!,
AiyagariProblem,
DynamicAiyagariProblem,
BansalYaronProblem,
BansalYaronProblemNewton,
BansalYaronProblemPowell

##############################################################################
##
## Load files
##
##############################################################################
include("aiyagari/aiyagari.jl")
include("aiyagari/dynamicaiyagari.jl")

include("bansalyaron/utils.jl")
include("bansalyaron/bansalyaronnewton.jl")
include("bansalyaron/bansalyaronpowell.jl")

end