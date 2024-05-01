using DiffGrating1D
using Plots; plotlyjs()
const DF = DiffGrating1D
include("geometry.jl")

geo = geometry();

##
plot(geo.domains[1],tangent=false,normal=false)
plot!(geo.domains[2],tangent=false,normal=false)