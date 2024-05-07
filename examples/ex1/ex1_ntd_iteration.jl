using DiffGrating1D
using Plots; plotlyjs()
using LinearAlgebra
const DF = DiffGrating1D
include("geometry.jl")

geo = geometry(;homogeneous=true);
#DF.check_geometry(geo;ϵtol=1e-4)

## NTD maps
domain = geo.domains[1]
k = domain.k
α0 = geo.α0
β = sqrt(Complex(k^2-α0^2))
@assert β == geo.βtop[geo.Jmax+1]
γ = exp(im*geo.α0*geo.L)

sol(x) = exp(im*(α0*x[1] - β*x[2]))
sol_grad(x) = DF.Point2D(im*α0,-im*β)*sol(x)
∂sol∂n(x,n) = DF.rdot(sol_grad(x),n)

## test Q and Y matrices
# initial conditions
Q = -geo.Gbottom*(-1)     # REVERSE SIGN!!
Y = I(size(Q,1))
γ = exp(im*geo.α0*geo.L)
# generate data
top_qpoints = (DF.qpoint(domain,i) for i in DF.topboundary_indices(domain))
bottom_qpoints = (DF.qpoint(domain,i) for i in DF.bottomboundary_indices(domain))

u_top = [sol(DF.point(q)) for q in top_qpoints]
∂u_top = [∂sol∂n(DF.point(q),DF.normal(q)) for q in top_qpoints]
u_bottom = [sol(DF.point(q)) for q in bottom_qpoints]
∂u_bottom = [∂sol∂n(DF.point(q),DF.normal(q)) for q in bottom_qpoints]

Qtop,Ytop = DF.obtain_Q_Y_matrices(domain,γ,Q,Y)

err_q1_inv = DF.rel_error(Qtop \ ∂u_top,u_top)  # works! 
err_y = DF.rel_error(Ytop*u_top,u_bottom)   # works!

## iterate to next domain
dom_2 = geo.domains[2]
top_top_qpoints = (DF.qpoint(dom_2,i) for i in DF.topboundary_indices(dom_2))
u_top_top = [sol(DF.point(q)) for q in top_top_qpoints]
∂u_top_top = [∂sol∂n(DF.point(q),DF.normal(q)) for q in top_top_qpoints]

Qtoptop_prev = -Qtop  # REVERSE SIGN!!
Qtoptop,Ytoptop = DF.obtain_Q_Y_matrices(dom_2,γ,Qtoptop_prev,Ytop)

err_q2_inv = DF.rel_error(Qtoptop \ ∂u_top_top,u_top_top) # works!
err_y2 = DF.rel_error(Ytoptop*u_top_top,u_bottom)    # works!

## solve top bulk
Qbulk = Qtoptop   # maintain sign, since normal still points up
qpoints = (DF.qpoint(dom_2,i) for i in DF.topboundary_indices(dom_2))
xpoints = (DF.point(q) for q in qpoints)
u_incident = [sol(x) for x in xpoints]
∂u_incident = [∂sol∂n(DF.point(q),DF.normal(q)) for q in qpoints]
rhs = [-2*im*β*u for u in u_incident]
Lmatrix = Qbulk-geo.Gtop
utop_apprx = Lmatrix\rhs

err_top = DF.rel_error(utop_apprx,u_incident)   # works!!!!
#err_top_forward = Lmatrix*u_incident-rhs

## obtain reflection and transmission coeff.
u_reflected = utop_apprx-u_incident
u_transmitted = Ytoptop*utop_apprx

# Fourier coefficients
r_coeff = geo.Ttop*u_reflected   # works!
t_coeff = geo.Tbottom*u_transmitted  # works!

## test 'solve_diffraction_problem'
u_reflected0,r_coeff0,u_transmitted0,t_coeff0 = DF.solve_diffraction_problem(geo);

# error vs by-hand
u_reflected0 ≈ u_reflected
r_coeff0 ≈ r_coeff
u_transmitted0 ≈ u_transmitted
t_coeff0 ≈ t_coeff

# error vs true
norm(u_reflected0,Inf)
DF.rel_error(u_transmitted0,u_bottom)

##
DF.check_homogeneous_geometry(geo;ϵtol=1e-4)
