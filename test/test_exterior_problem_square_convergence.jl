using DiffGrating1D
using LinearAlgebra
using Plots
const DF = DiffGrating1D

φ(t) = DF._curve_square(t)

Nlist = [40,80,120,160,200,240] .÷ 4
plist = [2,4,6,7,8]
k = 10
x0 = DF.Point2D(5.0,-1.0)     # exterior eval point

# exact solution
sol(x) = DF.hankel0(k*norm(x))
sol_x0 = sol(x0)

## test u = D[ϕd]
data = Dict()
err_noregularization = Float64[]
for N in Nlist
    for p in plist
        domain = DF.DomainSquare(k,N,p)
        sol_φ = [sol(DF.point(q)) for q in domain.quad]
        rhs = 2*sol_φ
        # dpot
        Zmatrix = DF.double_layer_matrix_plus_identity(domain)
        ϕ = Zmatrix\rhs
        Dpot = DF.double_layer_potential(domain,ϕ)
        sol_approx = Dpot(x0)
        error = abs((sol_x0-sol_approx)/sol_x0)
        data[(N,p)] = error
    end
    # compare with no regularization
    domain = DF.Domain(φ,k,4*N)  # times 4, bc there's only one curve
    sol_φ = [sol(DF.point(q)) for q in domain.quad]
    rhs = 2*sol_φ
    # dpot
    Zmatrix = DF.double_layer_matrix_plus_identity(domain)
    ϕ = Zmatrix\rhs
    Dpot = DF.double_layer_potential(domain,ϕ)
    sol_approx = Dpot(x0)
    error = abs((sol_x0-sol_approx)/sol_x0)
    push!(err_noregularization,error)
end

## Plots
fig = plot(xaxis=:log,yaxis=:log,legend=:bottomleft)
for p in plist
    errlist = [data[(N,p)] for N in Nlist]
    plot!(Nlist,errlist,label="p = $p")
end
plot!(Nlist,err_noregularization,label="no regularization",color=:black)
plot(fig)
