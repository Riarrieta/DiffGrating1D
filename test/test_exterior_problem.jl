using DiffGrating1D
using LinearAlgebra
using Plots
const DF = DiffGrating1D

φ(t) = DF.Point2D(cos(t)+0.65*cos(2*t)-0.65, 1.5*sin(t))

k = 1.5
N = 50
domain = DF.Domain(φ,k,N)
x0 = DF.Point2D(5.0,-1.0)     # exterior eval point

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

## test u = S[ϕ]
sol(x) = DF.hankel0(k*norm(x))
sol_x0 = sol(x0)
sol_φ = [sol(DF.point(q)) for q in domain.quad]

rhs = 2*sol_φ
Zmatrix = DF.single_layer_matrix(domain)
ϕ = Zmatrix\rhs

Spot = DF.single_layer_potential(domain,ϕ)
sol_approx = Spot(x0)
error = abs((sol_x0-sol_approx)/sol_x0)
