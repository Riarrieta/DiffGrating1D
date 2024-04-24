using DiffGrating1D
using LinearAlgebra
using ForwardDiff
using Plots
const DF = DiffGrating1D

φ(t) = DF.Point2D(cos(t), sin(t))

N = 50
k = 2.4048255576957727  # zero of J0
domain = DF.Domain(φ,k,N)

y0 = DF.Point2D(0.1,-0.1)     # exterior eval point
x0 = DF.Point2D(5,-1.0)     # exterior eval point
D(x,τ,k) = -0.5*DF.double_layer_kernel(x[1],x[2],τ,k)  # double layer kernel
S(x,τ,k) = 0.5*DF.single_layer_kernel(x[1],x[2],τ,k)  # single layer kernel

## interior test u = J0(k|x|), single layer
# exact solution
sol(x) = DF.besselj0(k*norm(x))
sol_grad(x) = -k*DF.besselj1(k*norm(x))*x/norm(x)
∂sol∂n(x,n) = dot(sol_grad(x),n)
sol_y0 = sol(y0)
#sol_φ = [sol(DF.point(q)) for q in domain.quad]
∂sol∂n_φ = [∂sol∂n(DF.point(q),DF.normal(q)) for q in domain.quad]

sol_approx = zero(ComplexF64)
for i in 1:N
    #u = sol_φ[i]
    ∂u∂n = ∂sol∂n_φ[i]
    q = domain.quad[i]
    sol_approx += S(y0,q,k)*∂u∂n
end
sol_approx = sol_approx*π/domain.n
error = abs((sol_y0-sol_approx)/sol_y0)

## interior test u = J0(k|x|), double layer
k2 = 3.83170597020751  # zero of J1
# exact solution
sol(x) = DF.besselj0(k2*norm(x))
sol_grad(x) = -k2*DF.besselj1(k2*norm(x))*x/norm(x)
∂sol∂n(x,n) = dot(sol_grad(x),n)
sol_y0 = sol(y0)
sol_φ = [sol(DF.point(q)) for q in domain.quad]
#∂sol∂n_φ = [∂sol∂n(DF.point(q),DF.normal(q)) for q in domain.quad]

sol_approx = zero(ComplexF64)
for i in 1:N
    u = sol_φ[i]
    #∂u∂n = ∂sol∂n_φ[i]
    q = domain.quad[i]
    sol_approx += -D(y0,q,k2)*u
end
sol_approx = sol_approx*π/domain.n
error = abs((sol_y0-sol_approx)/sol_y0)

## interior test u = planewave
# exact solution
#sol(x) = exp(im*k*x[1])
#sol_grad(x) = DF.Point2D(im*k*sol(x),0.0)
sol(x) = DF.besselj0(k*norm(x-x0))
sol_grad(x) = -k*DF.besselj1(k*norm(x-x0))*(x-x0)/norm(x-x0)
∂sol∂n(x,n) = dot(sol_grad(x),n)
sol_y0 = sol(y0)
sol_φ = [sol(DF.point(q)) for q in domain.quad]
∂sol∂n_φ = [∂sol∂n(DF.point(q),DF.normal(q)) for q in domain.quad]

sol_approx = zero(ComplexF64)
for i in 1:N
    u = sol_φ[i]
    ∂u∂n = ∂sol∂n_φ[i]
    q = domain.quad[i]
    sol_approx += S(y0,q,k)*∂u∂n - D(y0,q,k)*u
end
sol_approx = sol_approx*π/domain.n
error = abs((sol_y0-sol_approx)/sol_y0)