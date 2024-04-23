using DiffGrating1D
using ForwardDiff
using Plots
const DF = DiffGrating1D

φ(t) = DF.Point2D(cos(t)+0.65*cos(2*t)-0.65, 1.5*sin(t))

N = 20
tarray = range(0,2π-2π/N,N)

curve = φ.(tarray)
curve_x = [c[1] for c in curve]
curve_y = [c[2] for c in curve]

plot()
plot!(curve_x,curve_y,aspect_ratio=:equal)

##
a = DF._obtain_curve_parametrization(φ,tarray)
