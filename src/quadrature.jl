
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
    # derivative of change of variable
    # not used if change of variable is identity
    wprime::Float64
end
param(q::QPoint) = q.t
point(q::QPoint) = Point2D(q.x,q.y)
tvector(q::QPoint) = Point2D(q.tx,q.ty)
tnorm(q::QPoint) = q.tnorm
tnorm2(q::QPoint) = tnorm(q)^2
normal(q::QPoint) = Point2D(q.nx,q.ny)
ttvector(q::QPoint) = q.ttx,q.tty
distance(q1::QPoint,q2::QPoint) = sqrt((q1.x-q2.x)^2+(q1.y-q2.y)^2)
wprime(q::QPoint) = q.wprime

abstract type AbstractDomain end

struct Domain <: AbstractDomain
    n::Int64     # number of points = 2n
    k::Float64              # wavenumber
    quad::Vector{QPoint}    # quadrature
    tarray       # parametrization's parameter
    φ            # parametrization
end
wavenumber(d::Domain) = d.k
nunknowns(d::Domain) = 2*d.n
qpoint(d::Domain,i::Integer) = d.quad[i]

function Domain(φ,k,N)
    @assert iseven(N)
    n = N ÷ 2
    tarray = range(0,2π-2π/N,N)    # parameter in [0,2π)
    φp_func(t)  = ForwardDiff.derivative(φ,t)       # first derivative
    φpp_func(t) = ForwardDiff.derivative(φp_func,t) # second derivative
    # construct grid
    quad = QPoint[]
    for t in tarray
        x,y = φ(t)
        # tangent vector
        tx,ty = φp_func(t)
        tnorm = sqrt(tx^2+ty^2)
        # normal vector
        nx = ty/tnorm
        ny = -tx/tnorm
        # second derivative
        ttx,tty = φpp_func(t)
        wprime = NaN   # set to NaN, since it is not used here
        q = QPoint(t,x,y,tx,ty,tnorm,nx,ny,ttx,tty,wprime)
        push!(quad,q)
    end
    @assert length(quad) == N
    return Domain(n,k,quad,tarray,φ)
end

struct DomainWith1Corner <: AbstractDomain
    n::Int64     # number of points = 2n
    k::Float64              # wavenumber
    quad::Vector{QPoint}    # quadrature
    tarray       # parametrization's parameter
    φ            # parametrization
end
wavenumber(d::DomainWith1Corner) = d.k
nunknowns(d::DomainWith1Corner) = 2*d.n
qpoint(d::DomainWith1Corner,i::Integer) = d.quad[i]

# graded mesh change of variable
_vfunc(s,p) = (1/p-1/2)*((π-s)/π)^3 + 1/p*(s-π)/π+1/2
_wfunc(s,p) = 2π*_vfunc(s,p)^p/(_vfunc(s,p)^p + _vfunc(2π-s,p)^p)
_∂wfunc(s,p) = ForwardDiff.derivative((x) -> _wfunc(x,p),s)

function DomainWith1Corner(φ,k,N,p)
    @assert iseven(N)
    @assert p ≥ 2
    n = N ÷ 2
    tarray = range(0,2π-2π/N,N)    # parameter in [0,2π)
    φp_func(t)  = ForwardDiff.derivative(φ,t)       # first derivative
    φpp_func(t) = ForwardDiff.derivative(φp_func,t) # second derivative
    # construct grid
    quad = QPoint[]
    for t in tarray
        w = _wfunc(t,p)
        ∂w = _∂wfunc(t,p)
        # the first node derivative (the corner) is never used, set to NaN just in case
        if iszero(t)
            ∂w = NaN
        end
        x,y = φ(w)
        # tangent vector
        tx,ty = φp_func(w)
        tnorm = sqrt(tx^2+ty^2)
        # normal vector 
        nx = ty/tnorm
        ny = -tx/tnorm
        # second derivative
        ttx,tty = φpp_func(w)
        # we need to store the uniform t-parameter in the QPoint!
        q = QPoint(t,x,y,tx,ty,tnorm,nx,ny,ttx,tty,∂w)
        push!(quad,q)
    end
    @assert length(quad) == N
    return DomainWith1Corner(n,k,quad,tarray,φ)
end