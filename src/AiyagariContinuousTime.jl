module AiyagariContinuousTime

import StatsBase: denserank
import DataFrames: DataFrame
import Gadfly: plot, Geom
import Distances: chebyshev
using Gensys
using Calculus
using Gensys

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
simulate,
AiyagariProblem,
AiyagariArrays,
AiyagariSolution,
DynamicAiyagariSolution,
BansalYaronProblem,
plot_ll
#phact_solver

end