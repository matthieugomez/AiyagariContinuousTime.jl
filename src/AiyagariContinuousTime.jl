module AiyagariContinuousTime
import Distances: chebyshev
using Gensys
using Calculus

##############################################################################
##
## Load files
##
##############################################################################
include("aiyagari.jl")
include("dynamicaiyagari.jl")
##############################################################################
##
## Exported methods and types 
##
##############################################################################
export solve,
AiyagariProblem,
AiyagariArrays,
AiyagariSolution,
DynamicAiyagariSolution,
simulate
end