using Plots
using DiffGrating1D
using LinearAlgebra
const DF = DiffGrating1D

T0 = 1.0   # period
J0 = 7  # maximum mode
α0 = 0.3*π/T0 # Block wavenumber
max_freq = J0/T0
Jrange = -J0:J0
tvec = [min(1,1/abs(j^2))*rand(ComplexF64) for j in Jrange]
bvec = [min(1,1/abs(j^2))*rand(ComplexF64) for j in Jrange]

α(j) = α0+2π*j/T0

f(x) = sum(t*exp(im*α(j)*x) for (t,j) in zip(tvec,Jrange))

## plot
N = 15   # number of points
xrange = range(0,T0-T0/N,N)
fvec = f.(xrange)
#plot(xrange,@. real(fvec))

## error t
tjapprox = 1/T0*[sum(T0/N*fx*exp(-im*α(j)*x) for (x,fx) in zip(xrange,fvec)) for j in Jrange]
errt = DF.rel_error(tjapprox,tvec)

##  error g
g(x) = sum(t*b*exp(im*α(j)*x) for (t,b,j) in zip(tvec,bvec,Jrange))
gapprox(x) = sum(t*b*exp(im*α(j)*x) for (t,b,j) in zip(tjapprox,bvec,Jrange))

gvec = g.(xrange)
gapproxvec = gapprox.(xrange)
errg = DF.rel_error(gapproxvec,gvec)

## construct matrices
Tmatrix = [1/T0*T0/N*exp(-im*α(i)*xj) for i in Jrange, xj in xrange]
err_tmatrix = DF.rel_error(Tmatrix*fvec,tvec)
Fmatrix = [exp(im*α(j)*xi) for xi in xrange, j in Jrange]
err_fmatrix = DF.rel_error(Fmatrix*tvec,fvec)
Gmatrix = Fmatrix*diagm(bvec)*Tmatrix
err_gmatrix = DF.rel_error(Gmatrix*fvec,gvec)
Id_approx = Fmatrix*Tmatrix
err_identity = DF.rel_error(Id_approx*fvec,fvec)

## nonuniform mesh
N = 40   # number of points
p = 3
srange = range(0,2π-2π/N,N)  # uniform parameter in [0,2π]
wfunc(s) = T0/(2π)*DF._wfunc(s,p)
∂wfunc(s) = DF.ForwardDiff.derivative((x) -> wfunc(x),s)
xrange = wfunc.(srange)
∂wlist = ∂wfunc.(srange)
fvec = f.(xrange)

# error t
Tmatrix = [1/T0*2π/N*exp(-im*α(i)*xj)*∂wj for i in Jrange, (xj,∂wj) in zip(xrange,∂wlist)]
err_tmatrix = DF.rel_error(Tmatrix*fvec,tvec)

#error g
gvec = g.(xrange)
Fmatrix = [exp(im*α(j)*xi) for xi in xrange, j in Jrange]
Gmatrix = Fmatrix*diagm(bvec)*Tmatrix
err_gmatrix = DF.rel_error(Gmatrix*fvec,gvec)