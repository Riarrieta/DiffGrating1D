
hankel0(x) = hankelh1(0,x)  # Hankel function of first kind, order 0
hankel1(x) = hankelh1(1,x)  # Hankel function of first kind, order 1

# real dot product, without conjugation
rdot(a,b) = sum(i*j for (i,j) in zip(a,b))

const Point2D = SVector{2,T} where T
const Point3D = SVector{3,T} where T

function _obtain_curve_parametrization(φ,trange)
    φp_func(t)  = ForwardDiff.derivative(φ,t)
    φpp_func(t) = ForwardDiff.derivative(φp_func,t)

    φcurve = φ.(trange)
    φp = φp_func.(trange)
    normals = [Point2D(c[2],-c[1])/norm(c) for c in φp]
    φpp = φpp_func.(trange)
    return φcurve,φp,normals,φpp
end
