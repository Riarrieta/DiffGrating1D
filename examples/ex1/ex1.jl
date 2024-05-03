using DiffGrating1D
using Plots; plotlyjs()
using LinearAlgebra
const DF = DiffGrating1D
include("geometry.jl")

geo = geometry();
DF.check_geometry(geo;ϵtol=1e-4)
##
plot(geo.domains[1],tangent=false,normal=false)
plot!(geo.domains[2],tangent=false,normal=false)