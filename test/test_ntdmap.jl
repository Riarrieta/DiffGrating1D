using DiffGrating1D
using LinearAlgebra
using Plots
const DF = DiffGrating1D

φ(t) = DF.Point2D(cos(t)+0.65*cos(2*t)-0.65, 1.5*sin(t))

N = 50
Nhalf = N÷2
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

## use NtD map to solve for sol_φ[1:N/2] and ∂sol∂n_φ[1:N/2] 
# using sol_φ[N/2-1:end] and ∂sol∂n_φ[N/2-1:end]
D1 = @view Dmatrix[:,1:Nhalf]
D2 = @view Dmatrix[:,Nhalf+1:end]
S1 = @view Smatrix[:,1:Nhalf]
S2 = @view Smatrix[:,Nhalf+1:end]
sol_φ1 = @view sol_φ[1:Nhalf]
sol_φ2 = @view sol_φ[Nhalf+1:end]
∂sol∂n_φ1 = @view ∂sol∂n_φ[1:Nhalf]
∂sol∂n_φ2 = @view ∂sol∂n_φ[Nhalf+1:end]

Zmatrix = [D1 -S1]
rhs = -D2*sol_φ2 + S2*∂sol∂n_φ2
uapprox = Zmatrix \ rhs
utrue = vcat(sol_φ1,∂sol∂n_φ1)
sol_φ1_approx = uapprox[1:Nhalf]
∂sol∂n_φ1_approx = uapprox[Nhalf+1:end]
err1 = maximum(abs.(sol_φ1_approx-sol_φ1))/maximum(abs.(sol_φ1))
err2 = maximum(abs.(∂sol∂n_φ1_approx-∂sol∂n_φ1))/maximum(abs.(∂sol∂n_φ1))
err3 = maximum(abs.(Zmatrix*utrue-rhs))

## Dtn map
∂sol∂n_φ_approx = Smatrix \ (Dmatrix*sol_φ)
err_ntd = maximum(abs.(∂sol∂n_φ_approx-∂sol∂n_φ))/maximum(abs.(∂sol∂n_φ))

