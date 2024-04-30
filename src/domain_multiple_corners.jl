
struct DomainMultipleCorners <: AbstractDomainWithCorners
    n::Int64     # number of points = 2n
    k::Float64              # wavenumber
    quad::Vector{QPoint}    # quadrature
    corner_indices::Vector{Int64}  # global indices of corners
    edge_indices::Vector{Int64}    # global indices of edges
    tarray       # global parametrization's parameter
    φlist        # list of parametrization
end
wavenumber(d::DomainMultipleCorners) = d.k
nunknowns(d::DomainMultipleCorners) = 2*d.n
qpoint(d::DomainMultipleCorners,i::Integer) = d.quad[i]
corner_indices(d::DomainMultipleCorners) = d.corner_indices
corners(d::DomainMultipleCorners) = (qpoint(d,i) for i in corner_indices(d))
edge_indices(d::DomainMultipleCorners) = d.edge_indices
edges(d::DomainMultipleCorners) = (qpoint(d,i) for i in edge_indices(d))

# graded mesh change of variable
_ζfunc(s,si,sf) = (2*s-(si+sf))/(sf-si)
_w1func(s,p,si,sf) = (1/2-1/p)*(_ζfunc(s,si,sf))^3 + _ζfunc(s,si,sf)/p + 1/2
_w2func(s,p,si,sf) = 1-_w1func(s,p,si,sf)
_wgfunc(s,p,si,sf,ti,tf) = (tf*_w1func(s,p,si,sf)^p + ti*_w2func(s,p,si,sf)^p)/(_w1func(s,p,si,sf)^p+_w2func(s,p,si,sf)^p)
_∂wgfunc(p,si,sf,ti,tf) = (s) -> ForwardDiff.derivative((x)->_wgfunc(x,p,si,sf,ti,tf),s)

function _check_closed_curves(curves)
    # check if the collection of curves is closed
    # assume each curve is parametrized in [0,2π)
    check = true
    for i in 1:(length(curves)-1)
        φ1 = curves[i]
        φ2 = curves[i+1]
        check &= φ1(2π) ≈ φ2(0.0)
        !check && @info "" Tuple(φ1(2π)) Tuple(φ2(0.0)) norm(φ1(2π)-φ2(0.0))
    end
    # check last curve
    φ1 = curves[end]
    φ2 = curves[1]
    check &= φ1(2π) ≈ φ2(0.0)
    !check && @info "" Tuple(φ1(2π)) Tuple(φ2(0.0)) norm(φ1(2π)-φ2(0.0))
    return check
end
function DomainMultipleCorners(;φlist,k,Nlist,plist)
    @assert length(φlist) == length(Nlist) == length(plist)
    @assert _check_closed_curves(φlist)
    N = sum(Nlist)   # total number of nodes
    n = N ÷ 2
    tarray = range(0,2π-2π/N,N)    # global parameter in [0,2π)
    tindex = 1                     # gloabl parameter index
    # global indices of corners
    # first qnode of each curve is a corner
    corners_indices = vcat([1],Nlist[1:end-1])
    corners_indices = cumsum(corners_indices)
    # construct grid
    quad = QPoint[]
    #ti = 0.0
    #tf = 2π
    #si_index = 1
    #sf_index = 1
    for (φ,φN,p) in zip(φlist,Nlist,plist)
        @assert iseven(φN)
        @assert p ≥ 2
        φp_func(x)  = ForwardDiff.derivative(φ,x)       # first derivative
        φpp_func(x) = ForwardDiff.derivative(φp_func,x) # second derivative
        # update si,sf
        #si_index = sf_index
        #sf_index += φN  # sf_index belongs to the next curve
        # generate qpoints, assume φ is also parametrized in [0,2π)
        sarray = range(0,2π-2π/φN,φN)  # local parameter in [0,2π)
        for s in sarray
            w = _wfunc(s,p)
            ∂w = _∂wfunc(s,p)
            # the first node derivative (the corner) is never used, set to NaN just in case
            if iszero(s)
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
            # we need to store the global t-parameter in the QPoint!
            tglobal = tarray[tindex]
            tindex += 1
            q = QPoint(tglobal,x,y,tx,ty,tnorm,nx,ny,ttx,tty,∂w)
            push!(quad,q)
        end
        # update ti,tf
        #ti = tf
        #tf += 2π
    end
    @assert length(quad) == N
    @assert tindex == N+1
    edge_indices = setdiff(1:N,corners_indices)  # global indices of edges
    return DomainMultipleCorners(n,k,quad,corners_indices,edge_indices,tarray,φlist)
end

## Common shapes
function DomainSquare(k,N,p)
    φlist,Nlist,plist = curves_square(N,p)
    return DomainMultipleCorners(;φlist,k,Nlist,plist)
end
function DomainCosines(k,N,p)
    φlist,Nlist,plist = curves_cosines(N,p)
    return DomainMultipleCorners(;φlist,k,Nlist,plist)
end
