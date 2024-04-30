using DiffGrating1D
using LinearAlgebra
using Plots
const DF = DiffGrating1D

k = 1.5
N = 128
p = 8
domain = DF.DomainSquare(k,N,p)
#domain = DF.Domain(DF._curve_square,k,N)
#domain = DF.DomainCosines(k,N,p)
#domain = DF.Domain(DF._curve_cosines,k,N)
x0 = DF.Point2D(5.0,-1.0)     # exterior eval point

## test corner_modification_function(qi,j,d)
qj_index = 1
cm = [real(DF.corner_modification_function(domain.quad[i],qj_index,domain)) for i in 1:DF.nunknowns(domain)]
plot(cm)

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
