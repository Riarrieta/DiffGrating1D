using DiffGrating1D
using Plots; plotlyjs()
using LinearAlgebra
using JLD
using Interpolations
const DF = DiffGrating1D
include("geometry.jl")

# geometry
data_file = "examples/ex4/data_bic.jld"
Afrac = 0.1
Bfrac = 0.01
#data = Dict()   # Bfrac -> (wLlist,t0)
#@save data_file data
@load data_file data

## plot sample geometry
geoplot = geometry(;L=1,λ_over_L=1,Afrac,Bfrac,homogeneous=false);
plot(geoplot.domains[1],tangent=false,normal=false,aspect_ratio=:equal)
plot!(geoplot.domains[2],tangent=false,normal=false)

##
wLmin = 0.7719141 
wLmax = 0.7720703
current_wLlist = range(wLmin,wLmax,step=0.01/1024)
current_Blist = [0,Bfrac]
for B in current_Blist
    wLlist,t0list = get!(data,B,(Float64[],Float64[]))
    for wL in current_wLlist
        wL ∈ wLlist && continue
        @info "Computing wL = $wL, B = $B"
        λ_over_L = 1/wL   # wL/2πc = L/λ
        geo = geometry(;L=1,λ_over_L,Afrac,Bfrac=B,homogeneous=false);
        ## Diffraction problem
        _,_,_,t_coeff = DF.solve_diffraction_problem(geo);
        # diffraction efficiencies
        _,β0 = DF.αβfactors(geo)  # β0 of incident wave
        t_eff = @. geo.βbottom/β0*abs2(t_coeff)
        t0 = t_eff[geo.Jmax+1]  # transmission coeff zeroth order
        # add data
        push!(wLlist,wL)
        push!(t0list,t0)
        @save data_file data
    end
end

## Plot
@load data_file data
fig = plot(xlabel="w",ylabel="T₀")
data_keys = [0.0,Bfrac]
for B in data_keys
    wLlist,t0list = data[B]
    scatter!(wLlist,t0list,label="B = $B")
end
plot(fig)

## Plot exact transmittance for A = 0
@load data_file data
n1 = 1
n2 = sqrt(11.4)
n3 = 1
r12 = (n1-n2)/(n1+n2)
r23 = (n2-n3)/(n2+n3)
t12 = 2*n1/(n1+n2)
t23 = 2*n2/(n2+n3)
L = 1
d = 0.8*L/π
d0 = 0.1*L/π
D = d0+d/2  # thickness
ϕ(λ) = 2π*D*n2/λ
t(λ) = t12*t23*exp(-im*ϕ(λ))/(1+r12*r23*exp(-im*2*ϕ(λ)))
T(λ) = abs2(t(λ))

## Final plot
Tw(w) = T(1/w)
data_keys = [0.0,Bfrac]
fig = plot(xlabel="ωa/2πc",ylabel="T₀")
for B in data_keys
    wLlist,t0list = data[B]
    perm = sortperm(wLlist)
    wLlist = wLlist[perm]
    t0list = t0list[perm]
    # construct interpolant
    itp = interpolate((wLlist,), real.(t0list), Gridded(Linear()))
    plot!(wLlist,itp.(wLlist),label="A = $Afrac, B = $B")
end
# A = 0
wLlist0,_ = data[0.0]
wLlist0 = sort(wLlist0)
plot!(wLlist0,Tw.(wLlist0),label="A = 0.0, B = 0.0")
plot(fig)