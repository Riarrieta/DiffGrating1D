module DiffGrating1D

using LinearAlgebra
using StaticArrays
using SpecialFunctions
using ForwardDiff
using RecipesBase

include("utils.jl")
include("quadrature.jl")
include("domain_multiple_corners.jl")
include("kernels.jl")
include("matrix_discretization.jl")

end # module DiffGrating1D
