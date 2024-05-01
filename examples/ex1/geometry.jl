
function geometry(;L=1)
    ϵ1 = 1
    ϵ2 = 2.25
    p = 5
    λ = L/1.7
    k0 = 2π/λ
    θ = π/6
    d = L
    d0 = d/2
    d1 = d/2
    kB = k0*sqrt(ϵ1)
    kA = k0*sqrt(ϵ2)

    kvec1 = kB*DF.Point2D(sin(θ),-cos(θ))
    @assert kvec1 ≈ k0*DF.Point2D(1/2,-sqrt(3)/2)

    c0 = DF.Point2D(0.0,0.0)
    c1 = DF.Point2D(L,0.0)
    c2 = DF.Point2D(L,d0)
    c3 = DF.Point2D(0.0,d0)
    c4 = DF.Point2D(L,d0+d+d1)
    c5 = DF.Point2D(0.0,d0+d+d1)

    φ1 = DF.curve_straight_line(c0,c1)
    φ2 = DF.curve_straight_line(c1,c2)
    φ3(t) = DF.Point2D(L*(1-t/(2π)),-d/2*cos(2π*(1-t/(2π)))+(d0+d/2))
    φ4 = DF.curve_straight_line(c3,c0)
    φ6 = DF.curve_straight_line(c4,c2)
    φ7 = DF.curve_straight_line(c5,c4)
    φ8 = DF.curve_straight_line(c3,c5)
    N1 = 40
    N2 = 40
    N3 = 200
    N4 = 40
    N6 = 40
    N7 = 40
    N8 = 40

    # lower domain
    domA = DF.DomainMultipleCorners(;φlist=[φ1,φ2,φ3,φ4],
                                     k=kA,
                                     Nlist=[N1,N2,N3,N4],
                                     plist=[p,p,p,p],
                                     counterclockwise=true)
    # upper domain
    domB = DF.DomainMultipleCorners(;φlist=[φ3,φ8,φ7,φ6],
                                    k=kB,
                                    Nlist=[N3,N8,N7,N6],
                                    plist=[p,p,p,p],
                                    counterclockwise=false)
    # check common interface (φ3) has the same nodes
    domA_φ3_q = domA.quad[1+N1+N2 : N1+N2+N3]
    domB_φ3_q = domB.quad[1 : N3]
    @assert all(ai.x ≈ bi.x && ai.y ≈ bi.y for (ai,bi) in zip(domA_φ3_q,domB_φ3_q))

    return domA,domB
end