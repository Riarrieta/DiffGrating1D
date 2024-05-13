
function geometry(;nlevels=1,L=1,wL=0.35,θ=0.0,homogeneous=false)
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
    #θ = 0.0/360*2π  # angle in radians
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
