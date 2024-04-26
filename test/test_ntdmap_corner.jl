using DiffGrating1D
using LinearAlgebra
using Plots
const DF = DiffGrating1D

φ(t) = DF.curve_sharp_boomerang(t)

N = 128
p = 8
Nhalf = N÷2
k = 1.5
domain = DF.DomainWith1Corner(φ,k,N,p)
x0 = DF.Point2D(5.0,-1.0)     # exterior eval point

sol(x) = DF.hankel0(k*norm(x-x0))
sol_grad(x) = -k*DF.hankel1(k*norm(x-x0))*(x-x0)/norm(x-x0)
∂sol∂n(x,n) = DF.rdot(sol_grad(x),n)
sol_φ = [sol(DF.point(q)) for q in domain.quad]
∂sol∂n_φ = [∂sol∂n(DF.point(q),DF.normal(q)) for q in domain.quad]

## matrices
Dmatrix = DF.double_layer_matrix_plus_identity(domain)
Smatrix = DF.single_layer_matrix(domain)

## Ntd map
sol_φ_approx = Dmatrix \ (Smatrix*∂sol∂n_φ)
err_ntd = maximum(abs.(sol_φ_approx-sol_φ))/maximum(abs.(sol_φ))

## try Schur complement to eliminate corner variables
cv = [1]    # corner variables
ev = 2:N    # edge variables
D1 = Dmatrix[cv,cv]
D2 = Dmatrix[cv,ev]
D3 = Dmatrix[ev,cv]
D4 = Dmatrix[ev,ev]
S2 = Smatrix[cv,ev]
S4 = Smatrix[ev,ev]

D_D1 = (D4-D3*(D1\D2))  # Schur complement of D1
rhsmatrix = (S4-D3*(D1\S2))
rhs2 = rhsmatrix*∂sol∂n_φ[ev]
sol_φ_edge_approx = D_D1 \ rhs2
err_ntd_norcorner = maximum(abs.(sol_φ_edge_approx-sol_φ[ev]))/maximum(abs.(sol_φ[ev]))







