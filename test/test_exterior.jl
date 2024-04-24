using DiffGrating1D
using LinearAlgebra
using ForwardDiff
using Plots
const DF = DiffGrating1D

φ(t) = DF.Point2D(cos(t)+0.65*cos(2*t)-0.65, 1.5*sin(t))

N = 50
k = 1
domain = DF.Domain(φ,k,N)

# exact solution
x0 = DF.Point2D(5.0,-50.0)      # eval point
sol(x) = DF.hankel0(k*norm(x))
sol_grad(x) = -k*DF.hankel1(k*norm(x))*x/norm(x)
∂sol∂n(x,n) = dot(sol_grad(x),n)

sol_x0 = sol(x0)
sol_φ = [sol(DF.point(q)) for q in domain.quad]
∂sol∂n_φ = [∂sol∂n(DF.point(q),DF.normal(q)) for q in domain.quad]

## Test Green's representation exterior 1
D(x,τ,k) = -0.5*DF.double_layer_kernel(x[1],x[2],τ,k)  # double layer kernel
S(x,τ,k) = 0.5*DF.single_layer_kernel(x[1],x[2],τ,k)  # single layer kernel
sol_approx = zero(ComplexF64)
for i in 1:N
    u = sol_φ[i]
    ∂u∂n = ∂sol∂n_φ[i]
    q = domain.quad[i]
    sol_approx += D(x0,q,k)*u - S(x0,q,k)*∂u∂n
end
sol_approx = sol_approx*π/domain.n
error = abs((sol_x0-sol_approx)/sol_x0)

## 
# exact solution
x0 = DF.Point2D(5.0,0.0)      # eval point
sol(x) = exp(im*k*x[1])
sol_grad(x) = DF.Point2D(im*k*sol(x),0.0)
∂sol∂n(x,n) = dot(sol_grad(x),n)

sol_x0 = sol(x0)
sol_φ = [sol(DF.point(q)) for q in domain.quad]
∂sol∂n_φ = [∂sol∂n(DF.point(q),DF.normal(q)) for q in domain.quad]

# Test Green's representation exterior 2
D(x,τ,k) = -0.5*DF.double_layer_kernel(x[1],x[2],τ,k)  # double layer kernel
S(x,τ,k) = 0.5*DF.single_layer_kernel(x[1],x[2],τ,k)  # single layer kernel
sol_approx = zero(ComplexF64)
for i in 1:N
    u = sol_φ[i]
    ∂u∂n = ∂sol∂n_φ[i]
    q = domain.quad[i]
    sol_approx += D(x0,q,k)*u - S(x0,q,k)*∂u∂n
end
sol_approx = sol_approx*π/domain.n
error = abs((sol_x0-sol_approx)/sol_x0)

##
y0 = DF.Point2D([5.0,0.0])
sol2_re(x) = DF.besselj0(k*norm(x))
sol2_im(x) = DF.bessely0(k*norm(x))
sol2(x) = sol2_re(x) + im*sol2_im(x)
sol_grad2(x) = ForwardDiff.gradient(sol2_re,x) + im*ForwardDiff.gradient(sol2_im,x)
sol_grad2(y0)
sol_grad(y0)