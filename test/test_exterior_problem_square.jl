using DiffGrating1D
using LinearAlgebra
using Plots
const DF = DiffGrating1D

k = 1.5
N = 256
p = 8
domain = DF.DomainSquare(k,N,p)
x0 = DF.Point2D(5.0,-1.0)     # exterior eval point

## test kernels, N has to be large
q0index = 50
q0 = domain.quad[q0index]
q1 = domain.quad[q0index+1]
qnorm = sqrt((q0.x-q1.x)^2+(q0.y-q1.y)^2)
# M kernel
M1,M2 = DF.single_layer_kernel_M1_M2(q0,q0,k)
M1_approx,M2_approx = DF.single_layer_kernel_M1_M2(q1,q0,k)
# L kernel
L1,L2 = DF.double_layer_kernel_L1_L2(q0,q0,k)
L1_approx,L2_approx = DF.double_layer_kernel_L1_L2(q1,q0,k)

## test u = D[ϕ]
sol(x) = DF.hankel0(k*norm(x))
sol_x0 = sol(x0)
sol_φ = [sol(DF.point(q)) for q in domain.quad]

rhs = 2*sol_φ
Zmatrix = DF.double_layer_matrix_plus_identity(domain)
ϕ = Zmatrix\rhs

Dpot = DF.double_layer_potential(domain,ϕ)
sol_approx = Dpot(x0)
error = abs((sol_x0-sol_approx)/sol_x0)