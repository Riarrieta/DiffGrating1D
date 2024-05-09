
function geometry(;L=1,λ_over_L=1/1.9,Afrac=1,homogeneous=false)
    if !homogeneous
        # original structure
        ϵ1 = 1
        ϵ2 = 11.4
    else
        # set all ϵ to the maximum ϵ
        ϵ1 = 11.4
        ϵ2 = 11.4
    end
    p = 5
    λ = λ_over_L*L
    k0 = 2π/λ
    θ = 0/360*2π  # angle in radians
    d = 0.8*L/π
    d0 = 0.1*L/π
    d1 = d/2
    A = Afrac*d/2  # amplitude of sinusoid
    @assert d/2+d1>A && d/2+d0 >A
    ktop = k0*sqrt(ϵ1)   # air on the top
    ktopdom = ktop       # air on top domain
    kbottomdom = k0*sqrt(ϵ2)  # dielectric
    kbottom = ktop  # air on the bottom
    Jmax = 12

    kvec1 = ktop*DF.Point2D(sin(θ),-cos(θ))
    α0 = kvec1[1]

    c0 = DF.Point2D(0.0,0.0)
    c1 = DF.Point2D(L,0.0)
    c2 = DF.Point2D(L,d0)
    c3 = DF.Point2D(0.0,d0)
    c4 = DF.Point2D(L,d0+d+d1)
    c5 = DF.Point2D(0.0,d0+d+d1)

    φ1 = DF.curve_straight_line(c0,c1)
    φ2 = DF.curve_straight_line(c1,c2)
    φ3(t) = DF.Point2D(L*(1-t/(2π)),-A*cos(2π*(1-t/(2π)))+(d0+d/2))
    φ4 = DF.curve_straight_line(c3,c0)
    φ6 = DF.curve_straight_line(c4,c2)
    φ7 = DF.curve_straight_line(c5,c4)
    φ8 = DF.curve_straight_line(c3,c5)
    N1 = 120
    N2 = 60
    N3 = 120
    N4 = 60
    N6 = 60
    N7 = 120
    N8 = 60

    # lower domain
    domA_boundary_labels = [:bottom,:right,:top,:left]
    domA = DF.DomainMultipleCorners(;φlist=[φ1,φ2,φ3,φ4],
                                     k=kbottomdom,
                                     Nlist=[N1,N2,N3,N4],
                                     plist=[p,p,p,p],
                                     counterclockwise=true,
                                     φlabels=domA_boundary_labels)
    # upper domain
    domB_boundary_labels = [:bottom,:left,:top,:right]
    domB = DF.DomainMultipleCorners(;φlist=[φ3,φ8,φ7,φ6],
                                    k=ktopdom,
                                    Nlist=[N3,N8,N7,N6],
                                    plist=[p,p,p,p],
                                    counterclockwise=false,
                                    φlabels=domB_boundary_labels)

    geo = DF.Geometry([domA,domB],kbottom,ktop,L,α0,Jmax)
    return geo
end
