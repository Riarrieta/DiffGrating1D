using DiffGrating1D
using Plots; plotlyjs()
using LinearAlgebra
using JLD
const DF = DiffGrating1D
include("geometry.jl")

# First geometry
data_file = "examples/ex4/data1.jld"
#data = Dict()   # Afrac -> (λlist,t0)
#@save data_file data
@load data_file data

##
current_λlist = range(0.96875,1.0375,step=1/16)
current_Alist = [0.0,0.01,0.05]
for A in current_Alist
    λlist,t0list = get!(data,A,(Float64[],Float64[]))
    for λ in current_λlist
        λ ∈ λlist && continue
        @info "Computing λ = $λ, A = $A"
        geo = geo = geometry(;L=1,λ_over_L=λ,Afrac=A,homogeneous=false);
        ## Diffraction problem
        _,_,_,t_coeff = DF.solve_diffraction_problem(geo);
        # diffraction efficiencies
        _,β0 = DF.αβfactors(geo)  # β0 of incident wave
        t_eff = @. geo.βbottom/β0*abs2(t_coeff)
        t0 = t_eff[geo.Jmax+1]  # transmission coeff zeroth order
        # add data
        push!(λlist,λ)
        push!(t0list,t0)
        @save data_file data
    end
end

## Plot
@load data_file data
fig = plot(xlabel="λ/L",ylabel="T₀")
for A in keys(data)
    λlist,t0list = data[A]
    scatter!(λlist,t0list,label="A = $A")
end
plot(fig)

## Plot only one
@load data_file data
A_plot = 0.1
scatter(data[A_plot]...,label="A = $A_plot",xlabel="λ/L",ylabel="T₀")

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

A0 = 0.0
λlist0,t0list0 = data[A0]
λlist0_ord = range(0.5,1.5,200)
scatter(λlist0,t0list0,label="A = $A0",xlabel="λ/L",ylabel="T₀")
plot!(λlist0_ord,T.(λlist0_ord),label="true")