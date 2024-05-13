using DiffGrating1D
using Plots; plotlyjs()
using LinearAlgebra
using JLD
using Interpolations
const DF = DiffGrating1D
include("geometry.jl")

# Geometry
data_file = "examples/ex5/data_cond_number.jld"
wL = 1.25 # normalized freq wL/2πc
L = 1   # period
#clist = Float64[]
#αlist = Float64[]
#@save data_file αlist clist
@load data_file αlist clist

# Wood's anomalies
nmin = -ceil(2*wL)
nmax = ceil(wL)
αwlist = Float64[]
for n in nmin:nmax
    α1 = wL - n
    α2 = -wL - n
    (0 ≤ α1 ≤ wL) && α1 ∉ αwlist && push!(αwlist,α1)
    (0 ≤ α2 ≤ wL) && α2 ∉ αwlist && push!(αwlist,α2)
end
αwlist = sort(αwlist)

##
current_αlist = range(0.0, min(0.5,wL),step=2^-7)
for α in current_αlist
    α ∈ αlist && continue
    @info "Computing α = $α"
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
    push!(αlist,α)
    push!(clist,c)
    @info "Results:" α c size(Lmatrix)
    @save data_file data_file αlist clist
end

## Plot
@load data_file data_file αlist clist
scatter(αlist,clist,xlabel="",ylabel="")

## Pretty plot
@load data_file data_file αlist clist
perm = sortperm(αlist)
αlist_plot = αlist[perm]
clist_plot = clist[perm]
# construct interpolant
itp = interpolate((αlist_plot,), real.(clist_plot), Gridded(Linear()))
plot(αlist_plot,itp.(αlist_plot),xlabel="kₓa/2π",ylabel="cond. number",label="")
vline!(αwlist,label="Wood anomaly",legend=true,ls=:dash) 
xlims!(0,0.5)