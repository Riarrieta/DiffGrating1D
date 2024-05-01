using DiffGrating1D
using Plots; plotlyjs()
const DF = DiffGrating1D
include("geometry.jl")

domA,domB = geometry()

##
plot(domA,tangent=false,normal=false)
plot!(domB,tangent=false,normal=true)