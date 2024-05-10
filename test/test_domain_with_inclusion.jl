using DiffGrating1D
using LinearAlgebra
using ForwardDiff
using Plots
const DF = DiffGrating1D

bean(t) = 0.2*DF.curve_bean(t)

ext_N = 100
int_N = 100
w = 1
k = sqrt(2*w^2)
p = 5
ext_dom = DF.DomainSquare(k,ext_N,p)
int_dom = DF.Domain(bean,k,int_N)

y0 = DF.Point2D(0.5,0.5)     # interior eval point
x0 = DF.Point2D(0.0,0.0)     # interior interior eval point

sol(x) = cos(w*(x[1]))*cos(w*(x[2]))
sol_grad(x) = DF.ForwardDiff.gradient(sol,x)
∂sol∂n(x,n) = DF.rdot(sol_grad(x),n)
sol_y0 = sol(y0)
sol_x0 = sol(x0)

## plot
plot(ext_dom,tangent=false,normal=true,aspect_ratio=:equal)
plot!(int_dom,tangent=false,normal=true)

## test NtD exterior domain
uA = [sol(DF.point(q)) for q in ext_dom.quad]
∂uA = [∂sol∂n(DF.point(q),DF.normal(q)) for q in ext_dom.quad]

DpotA = DF.double_layer_potential(ext_dom,uA)
SpotA = DF.single_layer_potential(ext_dom,∂uA)
sol_approxA = SpotA(y0)-DpotA(y0)
errorA = abs((sol_y0-sol_approxA)/sol_y0)

DmatrixA = DF.double_layer_matrix_plus_identity(ext_dom)
SmatrixA = DF.single_layer_matrix(ext_dom)
uA_approx = DmatrixA \ (SmatrixA*∂uA)
errorA_ntdA = DF.rel_error(uA_approx,uA)

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
Dba = DF.double_layer_potential(int_dom,uB)
Sba = DF.single_layer_potential(int_dom,∂uB)
sol_approxAB = SpotA(y0) - DpotA(y0) - Sba(y0) + Dba(y0)
errorAB = abs((sol_y0-sol_approxAB)/sol_y0)

## test potential matrices
DAB = DF.double_layer_potential_matrix(ext_dom,int_dom)
SAB = DF.single_layer_potential_matrix(ext_dom,int_dom)
DBA = DF.double_layer_potential_matrix(int_dom,ext_dom)
SBA = DF.single_layer_potential_matrix(int_dom,ext_dom)

DAB_result = [-2*Daa(DF.point(q)) for q in int_dom.quad]  # add (-2) factor to correct kernel
SAB_result = [2*Saa(DF.point(q)) for q in int_dom.quad]   # add (2) factor to correct kernel
DBA_result = [-2*Dba(DF.point(q)) for q in ext_dom.quad]  # add (-2) factor to correct kernel
SBA_result = [2*Sba(DF.point(q)) for q in ext_dom.quad]   # add (2) factor to correct kernel

rab = 0.5*SAB_result-(-0.5)*DAB_result
DF.rel_error(rab,uB)

rba = 0.5*SBA_result-(-0.5)*DBA_result
norm(rba)

DAB_uA = DAB*uA
DAB_err = DF.rel_error(DAB_uA,DAB_result)

SAB_uA = SAB*∂uA
SAB_err = DF.rel_error(SAB_uA,SAB_result)

DBA_uB = DBA*uB
DBA_err = DF.rel_error(DBA_uB,DBA_result)

SBA_uB = SBA*∂uB
SBA_err = DF.rel_error(SBA_uB,SBA_result)

