
function geometry(;λ_over_L=1.7,Δd_over_L=0.4,L=1,homogeneous=false)
    if !homogeneous
        # original structure
        ϵ1 = 1
        ϵ2 = 5.0625
        ϵ3 = 2.1316
    else
        # set all ϵ to the maximum ϵ
        ϵ1 = 5.0625
        ϵ2 = 5.0625
        ϵ3 = 5.0625
    end
    p = 8
    λ = λ_over_L*L
    k0 = 2π/λ

    ktop = k0*sqrt(ϵ1)   # air on the top
    k1 = ktop       # air on top domain
    k2 = k0*sqrt(ϵ2)       # coating layer
    k3 = k0*sqrt(ϵ3)       # substrate
    kbottom = k3    # substrate on the bottom
    Jmax = 7
    θ = 0/360*2π  # angle in radians

    kvec1 = ktop*DF.Point2D(sin(θ),-cos(θ))
    α0 = kvec1[1]

    d = 2*L
    d0 = L/2
    d1 = d0
    Δd = Δd_over_L*L

    c0 = DF.Point2D(0.0,0.0)
    c1 = DF.Point2D(L,0.0)
    c2 = DF.Point2D(L,d0)
    c3 = DF.Point2D(L/2,d0+d)
    c4 = DF.Point2D(0.0,d0)

    c5 = c4 + DF.Point2D(0.0,Δd)
    c6 = c3 + DF.Point2D(0.0,Δd)
    c7 = c2 + DF.Point2D(0.0,Δd)

    c8 = c2 + DF.Point2D(0.0,d+Δd+d1)
    c9 = c4 + DF.Point2D(0.0,d+Δd+d1)

    φ1 = DF.curve_straight_line(c0,c1)
    φ2 = DF.curve_straight_line(c1,c2)
    φ3 = DF.curve_straight_line(c2,c3)
    φ4 = DF.curve_straight_line(c3,c4)
    φ5 = DF.curve_straight_line(c4,c0)

    φ6 = DF.curve_straight_line(c4,c5)
    φ7 = DF.curve_straight_line(c5,c6)
    φ8 = DF.curve_straight_line(c6,c7)
    φ9 = DF.curve_straight_line(c7,c2)

    φ10 = DF.curve_straight_line(c7,c8)
    φ11 = DF.curve_straight_line(c8,c9)
    φ12 = DF.curve_straight_line(c9,c5)

    N1 = 60
    N2 = 40
    N3 = 150
    N4 = 150
    N5 = 40

    N6 = 40
    N7 = 100
    N8 = 100
    N9 = 40

    N10 = 40
    N11 = 60
    N12 = 40

    # lower domain
    domA_boundary_labels = [:bottom,:right,:top,:top,:left]
    domA = DF.DomainMultipleCorners(;φlist=[φ1,φ2,φ3,φ4,φ5],
                                     k=k3,
                                     Nlist=[N1,N2,N3,N4,N5],
                                     plist=[p,p,p,p,p],
                                     counterclockwise=true,
                                     φlabels=domA_boundary_labels)
    # middle domain
    domB_boundary_labels = [:bottom,:bottom,:left,:top,:top,:right]
    domB = DF.DomainMultipleCorners(;φlist=[φ3,φ4,φ6,φ7,φ8,φ9],
                                    k=k2,
                                    Nlist=[N3,N4,N6,N7,N8,N9],
                                    plist=[p,p,p,p,p,p],
                                    counterclockwise=false,
                                    φlabels=domB_boundary_labels)
    # top domain
    domC_boundary_labels = [:bottom,:bottom,:right,:top,:left]
    domC = DF.DomainMultipleCorners(;φlist=[φ7,φ8,φ10,φ11,φ12],
                                    k=k1,
                                    Nlist=[N7,N8,N10,N11,N12],
                                    plist=[p,p,p,p,p],
                                    counterclockwise=true,
                                    φlabels=domC_boundary_labels)

    geo = DF.Geometry([domA,domB,domC],kbottom,ktop,L,α0,Jmax)
    return geo
end

function geometry_test(;λ_over_L=1.7,Δd_over_L=0.4,L=1,homogeneous=false)
    if !homogeneous
        # original structure
        ϵ1 = 1
        ϵ2 = 5.0625
        ϵ3 = 2.1316
    else
        # set all ϵ to the maximum ϵ
        ϵ1 = 5.0625
        ϵ2 = 5.0625
        ϵ3 = 5.0625
    end
    p = 8
    λ = λ_over_L*L
    k0 = 2π/λ

    ktop = k0*sqrt(ϵ1)   # air on the top
    k1 = ktop       # air on top domain
    k2 = k0*sqrt(ϵ2)       # coating layer
    k3 = k0*sqrt(ϵ3)       # substrate
    kbottom = k3    # substrate on the bottom
    Jmax = 7
    θ = 0/360*2π  # angle in radians

    kvec1 = ktop*DF.Point2D(sin(θ),-cos(θ))
    α0 = kvec1[1]

    d = 2*L
    d0 = L/2
    d1 = d0
    Δd = Δd_over_L*L

    c0 = DF.Point2D(0.0,0.0)
    c1 = DF.Point2D(L,0.0)
    c2 = DF.Point2D(L,d0)
    c3 = DF.Point2D(L/2,d0+d)
    c4 = DF.Point2D(0.0,d0)

    c5 = c4 + DF.Point2D(0.0,Δd)
    c6 = c3 + DF.Point2D(0.0,Δd)
    c7 = c2 + DF.Point2D(0.0,Δd)

    c8 = c2 + DF.Point2D(0.0,d+Δd+d1)
    c9 = c4 + DF.Point2D(0.0,d+Δd+d1)

    φ1 = DF.curve_straight_line(c0,c1)
    φ2 = DF.curve_straight_line(c1,c2)
    φ3 = DF.curve_straight_line(c2,c3)
    φ4 = DF.curve_straight_line(c3,c4)
    φ5 = DF.curve_straight_line(c4,c0)

    φ6 = DF.curve_straight_line(c4,c5)
    φ7 = DF.curve_straight_line(c5,c6)
    φ8 = DF.curve_straight_line(c6,c7)
    φ9 = DF.curve_straight_line(c7,c2)

    φ10 = DF.curve_straight_line(c7,c8)
    φ11 = DF.curve_straight_line(c8,c9)
    φ12 = DF.curve_straight_line(c9,c5)

    N1 = 40
    N2 = 40
    N3 = 100
    N4 = 100
    N5 = 40

    N6 = 40
    N7 = 100
    N8 = 100
    N9 = 40

    N10 = 40
    N11 = 40
    N12 = 40

    # lower domain
    domA_boundary_labels = [:bottom,:right,:top,:top,:left]
    domA = DF.DomainMultipleCorners(;φlist=[φ1,φ2,φ3,φ4,φ5],
                                     k=k3,
                                     Nlist=[N1,N2,N3,N4,N5],
                                     plist=[p,p,p,p,p],
                                     counterclockwise=true,
                                     φlabels=domA_boundary_labels)

    geo = DF.Geometry([domA],kbottom,ktop,L,α0,Jmax)
    return geo
end
