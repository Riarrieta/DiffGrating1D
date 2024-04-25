using DiffGrating1D
using ForwardDiff
using Plots
const DF = DiffGrating1D

#φ(t) = DF.Point2D(cos(t)+0.65*cos(2*t)-0.65, 1.5*sin(t))  # bean
φ(t) = DF.Point2D(2*sin(t/2), -sin(t))  # boomerang with corner

N = 50
domain = DF.Domain(φ,0.0,N)

curve_x = [q.x for q in domain.quad]
curve_y = [q.y for q in domain.quad]
φp_x = [q.tx for q in domain.quad]
φp_y = [q.ty for q in domain.quad]
normals_x = [q.nx for q in domain.quad]
normals_y = [q.ny for q in domain.quad]

## plotting
plot()
plot!(curve_x,curve_y,aspect_ratio=:equal)
quiver!(curve_x, curve_y, quiver=(φp_x, φp_y))
quiver!(curve_x, curve_y, quiver=(normals_x, normals_y))

## plot grade mesh change of variable
p = 7
nn = 50
sarray = range(0,2π-2π/nn,nn)
w = DF._wfunc.(sarray,p)
∂w = DF._∂wfunc.(sarray,p)
plot(sarray,w,label="w")
plot!(sarray,∂w,label="∂w")

## plot grade mesh
φ(t) = DF.Point2D(2*sin(t/2), -sin(t))  # boomerang with corner

N = 50
p = 5
domain = DF.DomainWith1Corner(φ,0.0,N,p)

curve_x = [q.x for q in domain.quad]
curve_y = [q.y for q in domain.quad]
φp_x = [q.tx for q in domain.quad]
φp_y = [q.ty for q in domain.quad]
normals_x = [q.nx for q in domain.quad]
normals_y = [q.ny for q in domain.quad]

## plotting
plot()
plot!(curve_x,curve_y,aspect_ratio=:equal)
quiver!(curve_x, curve_y, quiver=(φp_x, φp_y))
quiver!(curve_x, curve_y, quiver=(normals_x, normals_y))