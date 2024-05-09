using DiffGrating1D
using Plots; plotlyjs()
using LinearAlgebra
using JLD
const DF = DiffGrating1D
include("geometry.jl")

# First geometry
Δd_over_L = 0.4
#λlist = Float64[]
#rlist = Float64[]
#@save "examples/ex2/data.jld" λlist rlist
@load "examples/ex2/data.jld" λlist rlist

##
current_λlist = range(1.5375,1.5475,step=0.0025/2)
for λ in current_λlist
    λ ∈ λlist && continue
    @info "Computing λ = $λ"
    geo = geometry(;λ_over_L=λ,Δd_over_L,L=1,homogeneous=false);
    ## Diffraction problem
    _,r_coeff,_,_ = DF.solve_diffraction_problem(geo);
    # diffraction efficiencies
    _,β0 = DF.αβfactors(geo)  # β0 of incident wave
    r_eff = @. geo.βtop/β0*abs2(r_coeff)
    r0 = r_eff[geo.Jmax+1]  # reflection coeff zeroth order
    # add data
    push!(λlist,λ)
    push!(rlist,r0)
    @save "examples/ex2/data.jld" λlist rlist
end

## Plot
@load "examples/ex2/data.jld" λlist rlist
scatter(λlist,rlist,xlabel="λ/L",ylabel="R₀")