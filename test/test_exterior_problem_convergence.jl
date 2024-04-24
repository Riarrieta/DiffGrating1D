using DiffGrating1D
using LinearAlgebra
using Plots
const DF = DiffGrating1D

φ(t) = DF.Point2D(cos(t)+0.65*cos(2*t)-0.65, 1.5*sin(t))

Nlist = [40,80,120,160,200,320]
k = 20.0
x0 = DF.Point2D(5.0,-1.0)     # exterior eval point

# exact solution
sol(x) = DF.hankel0(k*norm(x))
sol_grad(x) = -k*DF.hankel1(k*norm(x))*x/norm(x)
sol_x0 = sol(x0)

## test u = D[ϕd] and u = S[ϕs]
edlist = Float64[]
eslist = Float64[]
for N in Nlist
    domain = DF.Domain(φ,k,N)

    sol_φ = [sol(DF.point(q)) for q in domain.quad]
    rhs = 2*sol_φ
    # dpot
    Zmatrix = DF.double_layer_matrix_plus_identity(domain)
    ϕ = Zmatrix\rhs
    Dpot = DF.double_layer_potential(domain,ϕ)
    sol_approx = Dpot(x0)
    error = abs((sol_x0-sol_approx)/sol_x0)
    push!(edlist,error)
    # spot
    Zmatrix = DF.single_layer_matrix(domain)
    ϕ = Zmatrix\rhs
    Spot = DF.single_layer_potential(domain,ϕ)
    sol_approx = Spot(x0)
    error = abs((sol_x0-sol_approx)/sol_x0)
    push!(eslist,error)
end

## Plots
plot(Nlist,eslist,xaxis=:log,yaxis=:log)
plot!(Nlist,edlist)
