using DiffGrating1D
using Plots; plotlyjs()
using LinearAlgebra
const DF = DiffGrating1D
include("geometry.jl")

geo = geometry();
#DF.check_geometry(geo;ϵtol=1e-4)

##
plot(geo.domains[1],tangent=false,normal=false)
plot!(geo.domains[2],tangent=false,normal=false)

## NTD maps
domain = geo.domains[1]
γ = exp(im*geo.α0*geo.L)
V = DF.obtain_ntd_map(domain)
N11,N12,N21,N22 = DF.obtain_reduced_ntd_map(domain,γ);

# test N matrices

## solve
# initial conditions
Q = -geo.Gbottom
Y = I(size(Q,1))
γ = exp(im*geo.α0*geo.L)
# iterate
for dom in geo.domains
    Q,Y = DF.obtain_Q_Y_matrices(dom,γ,Q,Y)
end
# solve top boundary
topdom = DF.topdomain(geo)
xpoints = (DF.qpoint(topdom,i).x for i in DF.topboundary_indices(topdom))
u_incident = [exp(im*geo.α0*x) for x in xpoints]
β0 = geo.βtop[geo.Jmax+1]
rhs = [-2*im*β0*u for u in u_incident]
Lmatrix = Q-geo.Gtop
utop = Lmatrix\rhs
u_reflected = utop-u_incident
u_transmitted = Y*utop

# Fourier coefficients
r_coeff = geo.Ttop*u_reflected
t_coeff = geo.Tbottom*u_transmitted

# diffraction efficiencies
r_eff = @. geo.βtop/β0*abs2(r_coeff)
t_eff = @. (geo.kbottom^2*geo.βbottom)/(geo.ktop^2*β0)*abs2(t_coeff)