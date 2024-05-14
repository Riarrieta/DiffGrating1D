using DiffGrating1D
using Plots; plotlyjs()
using LinearAlgebra
using JLD
using Interpolations
const DF = DiffGrating1D
include("geometry.jl")

# geometry
data_file = "examples/ex4/data_w_singlekx_T.jld"
Afrac = 0.5
dh = 0.0005 # step in both w̃ and α̃
α = 0.245
wLlist = range(0.5+dh,0.9,step=dh)
tlist = zeros(Float64,length(wLlist))
@save data_file α wLlist tlist
#@load data_file α wLlist tlist

##
progress = 1
total = length(wLlist)
Threads.@threads for i in 1:length(wLlist)
    wL = wLlist[i]
    @info "Computing α = $α, wL = $wL"        
    #@info "Progress $progress/$total"
    #progress +=1
    θ = asin(α/wL)
    λ_over_L = 1/wL    # wL/2πc = L/λ
    geo = geometry(;L=1,λ_over_L,Afrac,θ,homogeneous=false);
    ## Diffraction problem
    _,_,_,t_coeff = DF.solve_diffraction_problem(geo);
    # diffraction efficiencies
    _,β0 = DF.αβfactors(geo)  # β0 of incident wave
    t_eff = @. geo.βbottom/β0*abs2(t_coeff)
    t0 = t_eff[geo.Jmax+1]  # transmittance zeroth order
    # add data
    tlist[i] = t0
    #@save data_file α wLlist tlist
end
@save data_file α wLlist tlist

## Final plot
@load data_file α wLlist tlist
plot(wLlist,tlist)

