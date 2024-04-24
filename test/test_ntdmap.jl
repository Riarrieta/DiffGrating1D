using DiffGrating1D
using LinearAlgebra
using Plots
const DF = DiffGrating1D

φ(t) = DF.Point2D(cos(t)+0.65*cos(2*t)-0.65, 1.5*sin(t))

N = 50
k = 1.5
domain = DF.Domain(φ,k,N)

y0 = DF.Point2D(-0.1,0.3)      # interior eval point
x0 = DF.Point2D(5.0,-1.0)     # exterior eval point
D(x,τ,k) = -0.5*DF.double_layer_kernel(x[1],x[2],τ,k)  # double layer kernel
S(x,τ,k) = 0.5*DF.single_layer_kernel(x[1],x[2],τ,k)  # single layer kernel

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

## Dtn map
∂sol∂n_φ_approx = Smatrix \ (Dmatrix*sol_φ)
err_ntd = maximum(abs.(∂sol∂n_φ_approx-∂sol∂n_φ))/maximum(abs.(∂sol∂n_φ))

