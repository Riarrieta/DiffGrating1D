
struct QPoint
    # parameter t, such that φ(t) = (x,y), φ being the parametrization
    t::Float64
    # point
    x::Float64
    y::Float64
    # tangent vector
    tx::Float64
    ty::Float64
    tnorm::Float64
    # normal vector
    nx::Float64
    ny::Float64
    # second derivative
    ttx::Float64
    tty::Float64
end

param(q::QPoint) = q.t
point(q::QPoint) = q.x,q.y
tvector(q::QPoint) = q.tx,q.ty
tnorm(q::QPoint) = q.tnorm
tnorm2(q::QPoint) = tnorm(q)^2
nvector(q::QPoint) = q.nx,q.ny
ttvector(q::QPoint) = q.ttx,q.tty
distance(q1::QPoint,q2::QPoint) = sqrt((q1.x-q2.x)^2+(q1.y-q2.y)^2)

struct Domain
    k::Float64                 # wavenumber
    quad::Vector{QPoint}    # quadrature
end