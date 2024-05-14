using DiffGrating1D
using Plots; plotlyjs()
using LinearAlgebra
using JLD
using Interpolations
const DF = DiffGrating1D
include("geometry.jl")

# geometry
data_file = "examples/ex4/data_w_kx_T.jld"
Afrac = 0.5
dh = 0.005 # step in both w̃ and α̃
#αlist = range(0,0.5,step=dh)
#wLlist = range(0.5+dh,0.9,step=dh)
#tlist = zeros(Float64,length(αlist),length(wLlist))
#@save data_file αlist wLlist tlist
@load data_file αlist wLlist tlist

##
progress = 1
total = length(wLlist)*length(αlist)
Threads.@threads for j in 1:length(wLlist)
    wL = wLlist[j]
    for (i,α) in enumerate(αlist)
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
        tlist[i,j] = t0
        #@save data_file αlist wLlist tlist
    end
end
@save data_file αlist wLlist tlist

## Plot slice
@load data_file αlist wLlist tlist
index = 50
α = αlist[index]
tα = tlist[index,:]
plot(wLlist,tα)

## Final plot
@load data_file αlist wLlist tlist
heatmap(αlist,wLlist,transpose(tlist),xlabel="kₓa/2π",ylabel="ωa/2πc",c=:hot,clims=(0,1))
xlims!(0.0,0.5)
ylims!(0.5+dh,0.9)
