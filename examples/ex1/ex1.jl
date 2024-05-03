using DiffGrating1D
using Plots; plotlyjs()
using LinearAlgebra
const DF = DiffGrating1D
include("geometry.jl")

geo = geometry();

##
plot(geo.domains[1],tangent=false,normal=false)
plot!(geo.domains[2],tangent=false,normal=false)

## test boundary matrices bottom
domain = DF.bottomdomain(geo);boundary_label = :bottom;bvec = geo.βbottom
Tmatrix = geo.Tbottom; Gmatrix = geo.Gbottom
J0 = geo.Jmax
Jrange = -J0:J0
tvec = [min(1,1/abs(j^2))*rand(ComplexF64) for j in Jrange]
xvec = (DF.qpoint(domain,i).x for i in DF.boundary_indices(domain,boundary_label))
∂wlist = [DF.qpoint(domain,i).wprime for i in DF.boundary_indices(domain,boundary_label)]
npoints = length(xvec)
nfreq = length(Jrange)
f(x) = sum(tj*exp(im*αj*x) for (tj,αj) in zip(tvec,geo.αlist))
fvec = f.(xvec)
g(x) = sum(tj*im*bj*exp(im*αj*x) for (tj,bj,αj) in zip(tvec,bvec,geo.αlist))
gvec = g.(xvec)

tvec_approx = Tmatrix*fvec
dd = tvec_approx ./ tvec
@info "" DF.rel_error(tvec_approx,tvec)

gvec_approx = Gmatrix*fvec
dd2 = gvec_approx ./ gvec
@info "" DF.rel_error(gvec_approx,gvec)

## test boundary matrices top
domain = DF.topdomain(geo);boundary_label = :top; bvec = geo.βtop
Tmatrix = geo.Ttop; Gmatrix = geo.Gtop
J0 = geo.Jmax
Jrange = -J0:J0
tvec = [min(1,1/abs(j^2))*rand(ComplexF64) for j in Jrange]
xvec = (DF.qpoint(domain,i).x for i in DF.boundary_indices(domain,boundary_label))
npoints = length(xvec)
nfreq = length(Jrange)
f(x) = sum(tj*exp(im*αj*x) for (tj,αj) in zip(tvec,geo.αlist))
fvec = f.(xvec)
g(x) = sum(tj*im*bj*exp(im*αj*x) for (tj,bj,αj) in zip(tvec,bvec,geo.αlist))
gvec = g.(xvec)

tvec_approx = Tmatrix*fvec
dd = tvec_approx ./ tvec
@info "" DF.rel_error(tvec_approx,tvec)

gvec_approx = Gmatrix*fvec
dd2 = gvec_approx ./ gvec
@info "" DF.rel_error(gvec_approx,gvec)

## bottom, by hand
α(j) = geo.α0+2π*j/geo.L
alist0 = [α(j) for j in Jrange]
T0 = geo.L
N0 = length(DF.boundary_indices(domain,boundary_label))+1
φ1 = DF.curve_straight_line(DF.Point2D(0.0,0.0),DF.Point2D(T0,0.0))
srange = range(0,2π-2π/N0,N0)
p = 5
trange = DF._wfunc.(srange,p)
qvecs = φ1.(trange)
xvec0 = [q[1] for q in qvecs][2:end]
∂wlist0 = DF._∂wfunc.(srange,p)[2:end]
Tmatrix0 = [1/T0*2π/N0*exp(-im*α(i)*xj)*∂wj for i in Jrange, (xj,∂wj) in zip(xvec0,∂wlist0)]

tvec_approx2 = Tmatrix0/(2π)*fvec
dd0 = tvec_approx2 ./ tvec
@info "" DF.rel_error(tvec_approx2,tvec)
