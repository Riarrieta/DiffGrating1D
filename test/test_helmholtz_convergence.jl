using DiffGrating1D
using LinearAlgebra
using ForwardDiff
using Plots
const DF = DiffGrating1D

φ(t) = DF.Point2D(cos(t)+0.65*cos(2*t)-0.65, 1.5*sin(t))

Nlist = [30,40,50,60,70,80,90,100]
k = 10

x0 = DF.Point2D(5.0,-1.0)     # exterior eval point
D(x,τ,k) = -0.5*DF.double_layer_kernel(x[1],x[2],τ,k)  # double layer kernel
S(x,τ,k) = 0.5*DF.single_layer_kernel(x[1],x[2],τ,k)  # single layer kernel

# exact solution
sol(x) = DF.hankel0(k*norm(x))
sol_grad(x) = -k*DF.hankel1(k*norm(x))*x/norm(x)
∂sol∂n(x,n) = DF.rdot(sol_grad(x),n)
sol_x0 = sol(x0)

## Test Green's representation exterior, u = G
elist = Float64[]
for N in Nlist
    domain = DF.Domain(φ,k,N)

    sol_φ = [sol(DF.point(q)) for q in domain.quad]
    ∂sol∂n_φ = [∂sol∂n(DF.point(q),DF.normal(q)) for q in domain.quad]

    sol_approx = zero(ComplexF64)
    for i in 1:N
        u = sol_φ[i]
        ∂u∂n = ∂sol∂n_φ[i]
        q = domain.quad[i]
        sol_approx += D(x0,q,k)*u - S(x0,q,k)*∂u∂n
    end
    sol_approx = sol_approx*π/domain.n
    error = abs((sol_x0-sol_approx)/sol_x0)
    push!(elist,error)
end

## Plots
plot(Nlist, elist, xaxis=:log, yaxis=:log)



