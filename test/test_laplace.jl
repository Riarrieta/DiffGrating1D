using DiffGrating1D
using LinearAlgebra
using ForwardDiff
using Plots
const DF = DiffGrating1D

φ(t) = DF.Point2D(cos(t)+0.65*cos(2*t)-0.65, 1.5*sin(t))

N = 50
domain = DF.Domain(φ,0.0,N)

y0 = DF.Point2D(-0.1,0.3)      # interior eval point
x0 = DF.Point2D(5.0,-1.0)     # exterior eval point
D(x,τ) = DF.double_layer_kernel_laplace(x[1],x[2],τ)  # double layer kernel
S(x,τ) = DF.single_layer_kernel_laplace(x[1],x[2],τ)  # single layer kernel

## test laplace kernels
qindex = 12
q = domain.quad[qindex]
qtx,qty = DF.tvector(q)
qttx,qtty = DF.ttvector(q)
h = -0.001
qclose = DF.Point2D(q.x + h*qtx + h^2/2*qttx,q.y + h*qty + h^2/2*qtty)
H(x,τ) = DF.double_layer_kernel_laplace_correction(x[1],x[2],τ)  # double layer kernel for corner corrections
@assert H(x0,q) ≈ 2*D(x0,q)
@assert abs(H((q.x,q.y),q)-H(qclose,q)) < 1e-3
@assert abs(H((q.x,q.y),q)-2*D(qclose,q)) < 1e-3

## Interior test
# exact solution
sol(x) = 1.0
sol_y0 = sol(y0)
sol_φ = [sol(DF.point(q)) for q in domain.quad]

## Test Green's representation interior
sol_approx = zero(ComplexF64)
for i in 1:N
    u = sol_φ[i]
    q = domain.quad[i]
    sol_approx += -D(y0,q)*u
end
sol_approx = sol_approx*π/domain.n
error = abs((sol_y0-sol_approx)/sol_y0)

## Test Green's representation interior, for H kernel
qindex = 12
sol_q0 = -1   # should be -1 at the domain boundary
q0 = domain.quad[qindex]
sol_approx = zero(ComplexF64)
for i in 1:N
    q = domain.quad[i]
    sol_approx += H((q0.x,q0.y),q)*one(ComplexF64)
end
sol_approx = sol_approx*π/domain.n
error = abs((sol_q0-sol_approx)/sol_y0)

## exterior test u = 1
# exact solution
sol(x) = 1.0
sol_x0 = 0.0   # should be zero outside the domain
sol_φ = [sol(DF.point(q)) for q in domain.quad]

sol_approx = zero(ComplexF64)
for i in 1:N
    u = sol_φ[i]
    q = domain.quad[i]
    sol_approx += -D(x0,q)*u
end
sol_approx = sol_approx*π/domain.n
error = abs(sol_x0-sol_approx)

## exterior test u = G(x)
# exact solution
sol(x) = log(norm(x))
sol_grad(x) = 1/norm(x)^2*x
∂sol∂n(x,n) = rdot(sol_grad(x),n)
sol_x0 = sol(x0)
sol_φ = [sol(DF.point(q)) for q in domain.quad]
∂sol∂n_φ = [∂sol∂n(DF.point(q),DF.normal(q)) for q in domain.quad]

sol_approx = zero(ComplexF64)
for i in 1:N
    u = sol_φ[i]
    ∂u∂n = ∂sol∂n_φ[i]
    q = domain.quad[i]
    sol_approx += D(x0,q)*u - S(x0,q)*∂u∂n
end
sol_approx = sol_approx*π/domain.n
error = abs((sol_x0-sol_approx)/sol_x0)

## interior test u = G(x)
# exact solution
sol(x) = log(norm(x))
sol_grad(x) = 1/norm(x)^2*x
∂sol∂n(x,n) = rdot(sol_grad(x),n)
sol_y0 = 0.0   # should be zero inside the domain
sol_φ = [sol(DF.point(q)) for q in domain.quad]
∂sol∂n_φ = [∂sol∂n(DF.point(q),DF.normal(q)) for q in domain.quad]

sol_approx = zero(ComplexF64)
for i in 1:N
    u = sol_φ[i]
    ∂u∂n = ∂sol∂n_φ[i]
    q = domain.quad[i]
    sol_approx += D(y0,q)*u - S(y0,q)*∂u∂n
end
sol_approx = sol_approx*π/domain.n
error = abs(sol_y0-sol_approx)

