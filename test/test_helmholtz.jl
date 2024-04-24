using DiffGrating1D
using LinearAlgebra
using ForwardDiff
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

## test kernels
q = rand(domain.quad)
@assert D(x0,q,k) ≈ DF.double_layer_kernel_helmholtz(x0...,q,k)
@assert S(x0,q,k) ≈ DF.single_layer_kernel_helmholtz(x0...,q,k)

## test hankel functions
sol(x) = DF.hankel0(k*norm(x))
sol_grad(x) = -k*DF.hankel1(k*norm(x))*x/norm(x)
sol2_re(x) = DF.besselj0(k*norm(x))
sol2_im(x) = DF.bessely0(k*norm(x))
sol2(x) = sol2_re(x) + im*sol2_im(x)
sol_grad2(x) = ForwardDiff.gradient(sol2_re,x) + im*ForwardDiff.gradient(sol2_im,x)
@assert sol(x0) ≈ sol2(x0)
@assert sol_grad(x0) ≈ sol_grad2(x0)

# test Helmholtz equation
h = 0.0001
Δx = DF.Point2D(h,0)
Δy = DF.Point2D(0,h)
Δu = 1/h^2*(sol(x0+Δx)+sol(x0-Δx)+sol(x0+Δy)+sol(x0-Δy)-4*sol(x0))
Δu + k^2*sol(x0)

## Test Green's representation exterior, u = G
# exact solution
sol(x) = DF.hankel0(k*norm(x))
sol_grad(x) = -k*DF.hankel1(k*norm(x))*x/norm(x)
∂sol∂n(x,n) = DF.rdot(sol_grad(x),n)
sol_x0 = sol(x0)
sol_φ = [sol(DF.point(q)) for q in domain.quad]
∂sol∂n_φ = [∂sol∂n(DF.point(q),DF.normal(q)) for q in domain.quad]

Dpot = DF.double_layer_potential(domain,sol_φ)
Spot = DF.single_layer_potential(domain,∂sol∂n_φ)
sol_approx = Dpot(x0) - Spot(x0)
error = abs((sol_x0-sol_approx)/sol_x0)

## Test Green's representation exterior, u = planewave
# exact solution
sol(x) = exp(im*k*x[1])
sol_grad(x) = DF.Point2D(im*k*sol(x),0.0)
∂sol∂n(x,n) = DF.rdot(sol_grad(x),n)
sol_x0 = zero(ComplexF64)   # should be zero outside the domain
sol_φ = [sol(DF.point(q)) for q in domain.quad]
∂sol∂n_φ = [∂sol∂n(DF.point(q),DF.normal(q)) for q in domain.quad]

Dpot = DF.double_layer_potential(domain,sol_φ)
Spot = DF.single_layer_potential(domain,∂sol∂n_φ)
sol_approx = Dpot(x0) - Spot(x0)
error = abs((sol_x0-sol_approx)/sol_x0)
error = abs(sol_x0-sol_approx)

## Test Green's representation interior, u = G
# exact solution
sol(x) = DF.hankel0(k*norm(x-x0))
sol_grad(x) = -k*DF.hankel1(k*norm(x-x0))*(x-x0)/norm(x-x0)
∂sol∂n(x,n) = DF.rdot(sol_grad(x),n)
sol_y0 = sol(y0)
sol_φ = [sol(DF.point(q)) for q in domain.quad]
∂sol∂n_φ = [∂sol∂n(DF.point(q),DF.normal(q)) for q in domain.quad]

Dpot = DF.double_layer_potential(domain,sol_φ)
Spot = DF.single_layer_potential(domain,∂sol∂n_φ)
sol_approx = Spot(y0) - Dpot(y0)
error = abs((sol_y0-sol_approx)/sol_y0)
