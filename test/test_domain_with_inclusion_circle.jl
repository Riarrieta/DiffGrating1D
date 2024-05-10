using DiffGrating1D
using LinearAlgebra
using ForwardDiff
using Plots
const DF = DiffGrating1D

ext_N = 100
int_N = 100
kext = 2
kint = 2*kext
p = 5
r = 0.3 # circle radius
ext_dom = DF.DomainSquare(kext,ext_N,p;counterclockwise=false)
int_dom = DF.DomainCircle(kint,int_N,r)

y0 = DF.Point2D(0.5,0.5)     # interior eval point
x0 = DF.Point2D(0.0,0.0)     # interior interior eval point

# analytical solution
joext(x::Number) = DF.besselj0(kext*x)
joextp(x) = ForwardDiff.derivative(joext,x)
joint(x::Number) = DF.besselj0(kint*x)
jointp(x) = ForwardDiff.derivative(joint,x)

yoext(x::Number) = DF.bessely0(kext*x)
yoextp(x) = ForwardDiff.derivative(yoext,x)

brhs = [joint(r), jointp(r)]
Alhs = [joext(r) yoext(r);
        joextp(r) yoextp(r)]
A,B = Alhs \ brhs

sol(x) = A*DF.besselj0(kext*norm(x)) + B*DF.bessely0(kext*norm(x))
sol_grad(x) = DF.ForwardDiff.gradient(sol,x)
∂sol∂n(x,n) = DF.rdot(sol_grad(x),n)
sol_y0 = sol(y0)

# test sol
h = 0.001
Δx = DF.Point2D(h,0)
Δy = DF.Point2D(0,h)
eq = (sol(y0+Δx)+sol(y0-Δx)+sol(y0+Δy)+sol(y0-Δy)-4*sol(y0))/h^2+kext^2*sol(y0)

sol_int(x) = DF.besselj0(kint*norm(x))
sol_x0 = sol_int(x0)

uA = [sol(DF.point(q)) for q in ext_dom.quad]
∂uA = [∂sol∂n(DF.point(q),DF.normal(q)) for q in ext_dom.quad]

## plot
plot(ext_dom,tangent=false,normal=true,aspect_ratio=:equal)
plot!(int_dom,tangent=false,normal=true)

## test NtD exterior domain
DpotA = DF.double_layer_potential(ext_dom,uA)
SpotA = DF.single_layer_potential(ext_dom,∂uA)
sol_approxA = SpotA(y0)-DpotA(y0)
errorA = abs((sol_y0-sol_approxA)/sol_y0)   # it shouldn't work

DmatrixA = DF.double_layer_matrix_plus_identity(ext_dom)
SmatrixA = DF.single_layer_matrix(ext_dom)
uA_approx = DmatrixA \ (SmatrixA*∂uA)
errorA_ntdA = DF.rel_error(uA_approx,uA)   # it shouldn't work

## test NtD interior domain
uB = [sol(DF.point(q)) for q in int_dom.quad]
∂uB = [∂sol∂n(DF.point(q),DF.normal(q)) for q in int_dom.quad]

DpotB = DF.double_layer_potential(int_dom,uB)
SpotB = DF.single_layer_potential(int_dom,∂uB)
sol_approxB = SpotB(x0)-DpotB(x0)
errorB = abs((sol_x0-sol_approxB)/sol_x0)

DmatrixB = DF.double_layer_matrix_plus_identity(int_dom)
SmatrixB = DF.single_layer_matrix(int_dom)
uB_approx = DmatrixB \ (SmatrixB*∂uB)
errorB_ntdB = DF.rel_error(uB_approx,uB)

## test combined domain Green's formula
Daa = DF.double_layer_potential(ext_dom,uA)
Saa = DF.single_layer_potential(ext_dom,∂uA)
Dba = DF.double_layer_potential(int_dom,uB;k=ext_dom.k)
Sba = DF.single_layer_potential(int_dom,∂uB;k=ext_dom.k)
sol_approxAB = Saa(y0) - Daa(y0) - Sba(y0) + Dba(y0)
errorAB = abs((sol_y0-sol_approxAB)/sol_y0)

## potential matrices
DAB = DF.double_layer_potential_matrix(ext_dom,ext_dom,int_dom)
SAB = DF.single_layer_potential_matrix(ext_dom,ext_dom,int_dom)
DBA = DF.double_layer_potential_matrix(ext_dom,int_dom,ext_dom)
SBA = DF.single_layer_potential_matrix(ext_dom,int_dom,ext_dom)

##  check eq (1) and (2)
u1_approx = SmatrixA*∂uA-DmatrixA*uA-SBA*∂uB+DBA*uB
norm(u1_approx)

DmatrixB_kA = DF.double_layer_matrix_plus_identity(int_dom;k=ext_dom.k)
SmatrixB_kA = DF.single_layer_matrix(int_dom;k=ext_dom.k)
u2_approx = SAB*∂uA-DAB*uA-SmatrixB_kA*∂uB+DmatrixB_kA*uB
DF.rel_error(u2_approx,2*uB)

##  check eq (3)
N2 = DmatrixB \ SmatrixB  # NtD of interior domain
lhs_matrix = (2*N2+SmatrixB_kA-DmatrixB_kA*N2)
lhs_3 = lhs_matrix*∂uB
rhs_3 = SAB*∂uA-DAB*uA
DF.rel_error(lhs_3,rhs_3)

Z1 = lhs_matrix \ SAB
Z2 = lhs_matrix \ DAB
∂uB_approx_3 = Z1*∂uA-Z2*uA
err_∂u2 = DF.rel_error(∂uB,∂uB_approx_3)
NZ1 = N2*Z1
NZ2 = N2*Z2
uB_approx_3 = NZ1*∂uA-NZ2*uA
err_u2 = DF.rel_error(uB,uB_approx_3)

## check eq 4
lhs_matrix_4 = (DmatrixA - SBA*Z2 + DBA*NZ2)
rhs_matrix_4 = (SmatrixA - SBA*Z1 + DBA*NZ1)
lhs_4 = lhs_matrix_4*uA
rhs_4 = rhs_matrix_4*∂uA
DF.rel_error(lhs_4,rhs_4)
NtD2 = lhs_matrix_4 \ rhs_matrix_4

uA_approx_4 = lhs_matrix_4 \ rhs_4
DF.rel_error(uA_approx_4,uA)

norm(-SBA*Z2 + DBA*NZ2)
norm(-SBA*Z1 + DBA*NZ1)

## test obtain_pure_ntd_map(d::DomainWithInclusion)
dom = DF.DomainWithInclusion(ext_dom,int_dom)
N2 = DF.obtain_pure_ntd_map(dom)

uA_approx_5 = N2 * ∂uA
DF.rel_error(uA_approx_5,uA)
