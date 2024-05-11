using DiffGrating1D
using Plots; plotlyjs()
using LinearAlgebra
using JLD
using Interpolations
const DF = DiffGrating1D
include("geometry.jl")

# Geometry
data_file = "examples/ex5/data_multiple_wl.jld"
nlevels = 16
#wLlist = Float64[]
#tlist = Float64[]
#@save data_file wLlist tlist
@load data_file wLlist tlist

##
current_wLlist = range(0.4436524, 0.4453614,step=2^-14)
for wL in current_wLlist
    wL ∈ wLlist && continue
    @info "Computing wL = $wL"
    geo = geometry(;wL,nlevels=1,homogeneous=false);
    ## Diffraction problem
    _,_,_,t_coeff = DF.solve_diffraction_for_multiple_equal_domains(geo;nlevels);
    # diffraction efficiencies
    _,β0 = DF.αβfactors(geo)  # β0 of incident wave
    t_eff = @. geo.βbottom/β0*abs2(t_coeff)
    t0 = t_eff[geo.Jmax+1]  # transmittance zeroth order
    # add data
    push!(wLlist,wL)
    push!(tlist,t0)
    @info "Results:" wL t0
    @save data_file wLlist tlist
end

## Plot
@load data_file wLlist tlist
scatter(wLlist,tlist,xlabel="wL/2πc",ylabel="T₀")

## Pretty plot
@load data_file wLlist tlist
perm = sortperm(wLlist)
wLlist_plot = wLlist[perm]
t0list_plot = tlist[perm]
# construct interpolant
itp = interpolate((wLlist_plot,), real.(t0list_plot), Gridded(Linear()))
plot(wLlist_plot,itp.(wLlist_plot))