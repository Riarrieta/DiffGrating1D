using DiffGrating1D
using Plots; plotlyjs()
using LinearAlgebra
const DF = DiffGrating1D
include("geometry.jl")

# First geometry
geo = geometry(;λ_over_L=1.7,Δd_over_L=0.4,L=1,homogeneous=true);
DF._check_ntd(geo.domains[1],1,0,1,1e-4)
#DF.check_geometry(geo;ϵtol=1e-4)
#DF._check_T_G_matrices(geo,1e-4)
#DF._check_ntd(geo,1e-4)
#DF.check_homogeneous_geometry(geo;ϵtol=1e-4)

##
plot(geo.domains[1],tangent=false,normal=false,aspect_ratio=:equal)
plot!(geo.domains[2],tangent=false,normal=false)
plot!(geo.domains[3],tangent=false,normal=false)

## Diffraction problem
u_reflected,r_coeff,u_transmitted,t_coeff = DF.solve_diffraction_problem(geo);

# diffraction efficiencies
_,β0 = DF.αβfactors(geo)  # β0 of incident wave
r_eff = @. geo.βtop/β0*abs2(r_coeff)
t_eff = @. geo.βbottom/β0*abs2(t_coeff)   # (geo.kbottom^2*geo.βbottom)/(geo.ktop^2*β0)*abs2(t_coeff)

## Second geometry, L/λ = 1, θ = 30 degrees
geo2 = geometry_v2(;homogeneous=false);
#DF.check_geometry(geo2;ϵtol=1e-4)
#DF.check_homogeneous_geometry(geo2;ϵtol=1e-4)

# Diffraction problem
u_reflected2,r_coeff2,u_transmitted2,t_coeff2 = DF.solve_diffraction_problem(geo2);
# diffraction efficiencies
_,β02 = DF.αβfactors(geo2)  # β0 of incident wave
r_eff2 = @. geo2.βtop/β02*abs2(r_coeff2)
t_eff2 = @. geo2.βbottom/β02*abs2(t_coeff2)   # (geo.kbottom^2*geo.βbottom)/(geo.ktop^2*β0)*abs2(t_coeff)
