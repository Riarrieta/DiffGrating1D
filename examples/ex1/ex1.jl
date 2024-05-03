using DiffGrating1D
using Plots; plotlyjs()
using LinearAlgebra
const DF = DiffGrating1D
include("geometry.jl")

geo = geometry();
DF.check_geometry(geo;ϵtol=1e-4)
DF._check_ntd(geo,1e-4)
##
plot(geo.domains[1],tangent=false,normal=false)
plot!(geo.domains[2],tangent=false,normal=false)

##
geo_test = geometry_test(;L=1,NN=40,Ncurve=200);
DF._check_ntd(geo_test,1e-4)