using DiffGrating1D
using Plots; plotlyjs()
using LinearAlgebra
const DF = DiffGrating1D
include("geometry.jl")

# First geometry
wL = 0.5
nlevels = 36
geo = geometry(;wL,nlevels=1,homogeneous=false);
#DF.check_geometry(geo;ϵtol=1e-4)
#DF.check_homogeneous_geometry(geo;ϵtol=1e-4)

##
fig = plot(aspect_ratio=:equal)
for d in geo.domains
    plot!(d,tangent=false,normal=false)
end
plot(fig)

## Diffraction problem
u_reflected,r_coeff,u_transmitted,t_coeff = DF.solve_diffraction_for_multiple_equal_domains(geo;nlevels);


# diffraction efficiencies
_,β0 = DF.αβfactors(geo)  # β0 of incident wave
r_eff = @. geo.βtop/β0*abs2(r_coeff)
t_eff = @. geo.βbottom/β0*abs2(t_coeff)   # (geo.kbottom^2*geo.βbottom)/(geo.ktop^2*β0)*abs2(t_coeff)
t_eff[geo.Jmax+1],r_eff[geo.Jmax+1]  # transmitance/reflectance zeroth order