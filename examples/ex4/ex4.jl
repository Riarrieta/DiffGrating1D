using DiffGrating1D
using Plots; plotlyjs()
using LinearAlgebra
const DF = DiffGrating1D
include("geometry.jl")

# First geometry
geo = geometry(;L=1,λ_over_L=1/1.9,Afrac=0.9,homogeneous=false);
#DF.check_geometry(geo;ϵtol=1e-4)
#DF.check_homogeneous_geometry(geo;ϵtol=1e-4)

##
plot(geo.domains[1],tangent=false,normal=false,aspect_ratio=:equal)
plot!(geo.domains[2],tangent=false,normal=false)

## Diffraction problem
u_reflected,r_coeff,u_transmitted,t_coeff = DF.solve_diffraction_problem(geo);

# diffraction efficiencies
_,β0 = DF.αβfactors(geo)  # β0 of incident wave
r_eff = @. geo.βtop/β0*abs2(r_coeff)
t_eff = @. geo.βbottom/β0*abs2(t_coeff)   # (geo.kbottom^2*geo.βbottom)/(geo.ktop^2*β0)*abs2(t_coeff)
