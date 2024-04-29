
hankel0(x) = hankelh1(0,x)  # Hankel function of first kind, order 0
hankel1(x) = hankelh1(1,x)  # Hankel function of first kind, order 1

# real dot product, without conjugation
rdot(a,b) = sum(i*j for (i,j) in zip(a,b))

const Point2D = SVector{2,T} where T
const Point3D = SVector{3,T} where T

# curves
curve_bean(t) = Point2D(cos(t)+0.65*cos(2*t)-0.65, 1.5*sin(t))
curve_sharp_boomerang(t) = Point2D(2*sin(t/2)-1, -sin(t))  # boomerang with corner
curve_straight_line(p1::Point2D,p2::Point2D) = t -> p1 + t/(2π)*(p2-p1)
function curves_square(N,p)
    # square with corners (1,1), (-1,1), (-1,-1),(1,-1)
    @assert iseven(N)
    Nlist = [N for _ in 1:4]
    plist = [p for _ in 1:4]
    p1 = Point2D(1.0,1.0)
    p2 = Point2D(-1.0,1.0)
    p3 = Point2D(-1.0,-1.0)
    p4 = Point2D(1.0,-1.0)
    φ1 = curve_straight_line(p1,p2)
    φ2 = curve_straight_line(p2,p3)
    φ3 = curve_straight_line(p3,p4)
    φ4 = curve_straight_line(p4,p1)
    φlist = [φ1,φ2,φ3,φ4]
    return φlist,Nlist,plist
end

function _obtain_curve_parametrization(φ,trange)
    φp_func(t)  = ForwardDiff.derivative(φ,t)
    φpp_func(t) = ForwardDiff.derivative(φp_func,t)

    φcurve = φ.(trange)
    φp = φp_func.(trange)
    normals = [Point2D(c[2],-c[1])/norm(c) for c in φp]
    φpp = φpp_func.(trange)
    return φcurve,φp,normals,φpp
end
