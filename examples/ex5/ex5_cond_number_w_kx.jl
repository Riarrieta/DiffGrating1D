using DiffGrating1D
using Plots; plotlyjs()
using LinearAlgebra
using JLD
using Interpolations
const DF = DiffGrating1D
include("geometry.jl")

# Geometry
data_file = "examples/ex5/data_cond_number_w_kx.jld"
L = 1   # period
dh = 0.005 # step in both w̃ and α̃
#αlist = range(0,0.5,step=dh)
#wLlist = range(0.5,1.5,step=dh)
#clist = zeros(Float64,length(αlist),length(wLlist))
#@save data_file αlist wLlist clist
@load data_file αlist wLlist clist

##
Threads.@threads for j in 1:length(wLlist)
    wL = wLlist[j]
    for (i,α) in enumerate(αlist)
        if α ≥ wL
            clist[i,j] = NaN
            continue
        end
        @info "Computing α = $α, wL = $wL"  
        θ = asin(α/wL)
        geo = geometry(;L,wL,θ,nlevels=1,homogeneous=false);
        ## Compute system matrix
        # initial conditions
        γ = DF.γfactor(geo)
        Q = -geo.Gbottom  # -iB^(2)
        Y = I(size(Q,1))
        # iterate
        for (i,dom) in enumerate(geo.domains) # from bottom to top
            @info "Solving Domain $i"
            # multiply by -1, since the normals of adyacent domains point in opposite directions
            rmul!(Q,-one(ComplexF64))  
            # obtain Q and Y from Domain
            Q,Y = DF.obtain_Q_Y_matrices(dom,γ,Q,Y)
        end
        Lmatrix = Q-geo.Gtop  # Q-iB^(1)
        c = cond(Lmatrix)
        # add data
        clist[i,j] = c
        @info "Results:" c size(Lmatrix)
    end
end
@save data_file αlist wLlist clist

## Plot slice
@load data_file αlist wLlist clist
index = 50
α = αlist[index]
tα = clist[index,:]
plot(wLlist,tα)

## Final plot
@load data_file αlist wLlist clist
fig = heatmap(αlist,wLlist,(transpose(clist)),grid=false,xlabel="kₓa/2π",ylabel="ωa/2πc",c=:hot)
xlims!(0.0,0.5)
ylims!(0.5,1.5)
# +-(a_n + n) =w
#= plot Wood's anomalies
for n in -10:10
    plot!(αlist,αlist .+ n)
    plot!(αlist,-(αlist .+ n))
end
plot(fig) =#
