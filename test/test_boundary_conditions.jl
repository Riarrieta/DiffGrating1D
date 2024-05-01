using Plots
using DiffGrating1D

T0 = π   # period
J0 = 15  # maximum mode
α0 = 0.3*π/T0 # Block wavenumber
max_freq = J0/T0
Jrange = -J0:J0
tvec = [min(1,1/abs(j^2))*rand(ComplexF64) for j in Jrange]
bvec = [min(1,1/abs(j^2))*rand(ComplexF64) for j in Jrange]

α(j) = α0+2π*j/T0

f(x) = sum(t*exp(im*α(j)*x) for (t,j) in zip(tvec,Jrange))

## plot
N = 31   # number of points
xrange = range(0,T0-T0/N,N)
fvec = f.(xrange)
#plot(xrange,@. real(fvec))

##
tjapprox = 1/T0*[sum(T0/N*fx*exp(-im*α(j)*x) for (x,fx) in zip(xrange,fvec)) for j in Jrange]
errt = DF.rel_error(tjapprox,tvec)

##
g(x) = sum(t*b*exp(im*α(j)*x) for (t,b,j) in zip(tvec,bvec,Jrange))
gapprox(x) = sum(t*b*exp(im*α(j)*x) for (t,b,j) in zip(tjapprox,bvec,Jrange))

gvec = g.(xrange)
gapproxvec = gapprox.(xrange)
errg = DF.rel_error(gapproxvec,gvec)