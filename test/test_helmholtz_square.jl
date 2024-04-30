using DiffGrating1D
using LinearAlgebra
using ForwardDiff
using Plots
const DF = DiffGrating1D

φ(t) = DF._curve_square(t)
#φ(t) = DF.Point2D(cos(t)+0.65*cos(2*t)-0.65, 1.5*sin(t))

N = 1024
w = π/2
k = sqrt(2*w^2)
p = 8
#domain = DF.Domain(φ,k,N)
domain = DF.DomainSquare(k,N,p)

y0 = DF.Point2D(0.1,-0.1)     # exterior eval point
x0 = DF.Point2D(5.0,-1.0)     # exterior eval point

## test laplacian
u(x) = cos(w*x[1])*cos(w*x[2])
h = 0.01
Δx = DF.Point2D(h,0)
Δy = DF.Point2D(0,h)
eq = (u(y0+Δx)+u(y0-Δx)+u(y0+Δy)+u(y0-Δy)-4*u(y0))/h^2+k^2*u(y0)

## interior test u = cos(kx)*cos(ky)
# exact solution
#sol(x) = DF.besselj0(k*norm(x-x0))#cos(k*x[1])
sol(x) = cos(w*(x[1]-1))*cos(w*(x[2]-1))
#sol_grad(x) = DF.Point2D(-w*sin(w*x[1])*cos(w*x[2]),-w*cos(w*x[1])*sin(w*x[2]))
sol_grad(x) = DF.ForwardDiff.gradient(sol,x)
∂sol∂n(x,n) = DF.rdot(sol_grad(x),n)
sol_y0 = sol(y0)
sol_φ = [sol(DF.point(q)) for q in domain.quad]
∂sol∂n_φ = [∂sol∂n(DF.point(q),DF.normal(q)) for q in domain.quad]


Dpot = DF.double_layer_potential(domain,sol_φ)
Spot = DF.single_layer_potential(domain,∂sol∂n_φ)
sol_approx = -Dpot(y0)*4
error = abs((sol_y0-sol_approx)/sol_y0)

## interior test u = sin(kx)*sin(ky)
# exact solution
sol(x) = sin(w*(x[1]-1))*sin(w*(x[2]-1))
sol_grad(x) = DF.ForwardDiff.gradient(sol,x)
∂sol∂n(x,n) = DF.rdot(sol_grad(x),n)
sol_y0 = sol(y0)
sol_φ = [sol(DF.point(q)) for q in domain.quad]
∂sol∂n_φ = [∂sol∂n(DF.point(q),DF.normal(q)) for q in domain.quad]


Dpot = DF.double_layer_potential(domain,sol_φ)
Spot = DF.single_layer_potential(domain,∂sol∂n_φ)
sol_approx = Spot(y0)*4-Dpot(y0)*4
error = abs((sol_y0-sol_approx)/sol_y0)
