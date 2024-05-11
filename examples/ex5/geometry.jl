
function geometry(;nlevels=1,L=1,wL=0.35,homogeneous=false)
    if !homogeneous
        # original structure
        ϵ1 = 1
        ϵ2 = 8.9
    else
        # set all ϵ to the maximum ϵ
        ϵ1 = 8.9
        ϵ2 = 8.9
    end
    p = 5
    λ_over_L = 1/wL  #(wL/2πc = L/λ)
    λ = λ_over_L*L
    k0 = 2π/λ
    r = 0.2*L  # radius circle
    θ = 0.0/360*2π  # angle in radians
    Jmax = 7
    
    ktop = k0*sqrt(ϵ1)   # air on the top
    kext = ktop       # air in interior domain
    kint = k0*sqrt(ϵ2)  # dielectric
    kbottom = ktop  # air on the bottom

    kvec1 = ktop*DF.Point2D(sin(θ),-cos(θ))
    α0 = kvec1[1]

    Next = 50
    Nint = 100

    function domain(;level=0,counterclockwise)
        z0 = L/2
        z = z0 + L*level
        #println((level,z,counterclockwise))
        # exterior domain
        ext_dom = DF.DomainSquare(kext,Next,p;counterclockwise,L,z)
        # interior domain = 
        int_dom = DF.DomainCircle(kint,Nint,r;z)
        return DF.DomainWithInclusion(ext_dom,int_dom)
    end
    cc_func(level) = (-1)^level == 1
    domains = [domain(;level=i,counterclockwise=cc_func(i)) for i in 0:nlevels-1]
    geo = DF.Geometry(domains,kbottom,ktop,L,α0,Jmax)
    return geo
end

function solve_for_multiple_same_domain(geo::DF.Geometry;nlevels)
    @assert length(geo.domains) == 1
    @assert geo.domains[1] isa DF.DomainWithInclusion
    domain = geo.domains[1]
    # initial conditions
    γ = DF.γfactor(geo)
    α0,β0 = DF.αβfactors(geo)
    Q = -geo.Gbottom  # -iB^(2)
    Y = I(size(Q,1))
    # domain N matrices
    N11,N12,N21,N22 = DF.obtain_reduced_ntd_map(domain,γ)
    # iterate
    for (i) in 1:nlevels # from bottom to top
        @info "Solving Domain $i"
        # multiply by -1, since the normals of adyacent domains point in opposite directions
        rmul!(Q,-one(ComplexF64))  
        # obtain Q and Y from Domain
        Z = (I-N11*Q) \ N12
        Q = inv(N22+N21*Q*Z)
        Y = Y*Z*Q
    end
    @info "Solving Top Domain"
    # solve top boundary
    topdom = DF.topdomain(geo)
    top_points = (DF.point(DF.qpoint(topdom,i)) for i in DF.topboundary_indices(topdom))
    L = geo.L
    Δz = (nlevels-1)*L
    u_inc(x) = exp(im*(α0*x[1]-β0*(x[2]+Δz)))   # incident planewave
    u_incident = [u_inc(x) for x in top_points]   
    β0 = geo.βtop[geo.Jmax+1]      # β of incident planewave
    rhs = [-2*im*β0*u for u in u_incident]
    Lmatrix = Q-geo.Gtop  # Q-iB^(1)
    # solve and get fields
    utop = Lmatrix\rhs
    u_reflected = utop-u_incident
    u_transmitted = Y*utop
    # get reflection and transmission coeff
    r_coeff = geo.Ttop*u_reflected
    t_coeff = geo.Tbottom*u_transmitted
    @info "Done"
    return u_reflected,r_coeff,u_transmitted,t_coeff
end