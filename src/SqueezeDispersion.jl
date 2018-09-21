module SqueezeDispersion

using FourierFlows, PassiveTracerFlows, PyPlotPlus
using FFTW, Printf, PyPlot, Parameters, JLD2

import PassiveTracerFlows.TracerAdvDiff
import LinearAlgebra: mul!, ldiv!

export setupproblem, runproblem, rundefaultproblem, sinusoidalbathymetrymask, constantkappaenhancement
include("undulatingbarotropicflow.jl")

export genfig2data, plotfig2enhancement
include("manuscripttools.jl")

export squeezediffusivity, wavychannelmask, getkappa, getgrid, getgridspacing, getgrid, getxy,  calcenhancement
include("analysis.jl")

# Samoan module
include("samoan.jl")

end # module
