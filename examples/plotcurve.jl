using DiffGrating1D
using ForwardDiff
using Plots
const DF = DiffGrating1D

φ(t) = DF.Point2D(cos(t)+0.65*cos(2*t)-0.65, 1.5*sin(t))

N = 50
tarray = range(0,2π-2π/N,N)

φcurve,φp,normals,φpp = DF._obtain_curve_parametrization(φ,tarray)
curve_x = [c[1] for c in φcurve]
curve_y = [c[2] for c in φcurve]
φp_x = [c[1] for c in φp]
φp_y = [c[2] for c in φp]
normals_x = [c[1] for c in normals]
normals_y = [c[2] for c in normals]

## plotting

plot()
plot!(curve_x,curve_y,aspect_ratio=:equal)
quiver!(curve_x, curve_y, quiver=(φp_x, φp_y))
quiver!(curve_x, curve_y, quiver=(normals_x, normals_y))
