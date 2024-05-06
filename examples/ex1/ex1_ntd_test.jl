using DiffGrating1D
using Plots; plotlyjs()
using LinearAlgebra
const DF = DiffGrating1D
include("geometry.jl")

geo = _simple_geometry();
#DF.check_geometry(geo;ϵtol=1e-4)

## NTD maps
domain = geo.domains[1]
k = domain.k
α0 = geo.α0
β = sqrt(Complex(k^2-α0^2))
γ = exp(im*geo.α0*geo.L)

sol(x) = exp(im*(α0*x[1] - β*x[2]))
sol_grad(x) = DF.Point2D(im*α0,-im*β)*sol(x)
∂sol∂n(x,n) = DF.rdot(sol_grad(x),n)
sol_φ = [sol(DF.point(q)) for q in domain.quad]
∂sol∂n_φ = [∂sol∂n(DF.point(q),DF.normal(q)) for q in domain.quad]

## test full NtD map
Dmatrix = DF.double_layer_matrix_plus_identity(domain)
Smatrix = DF.single_layer_matrix(domain)
# Ntd map
sol_φ_approx = Dmatrix \ (Smatrix*∂sol∂n_φ)
err_ntd = DF.rel_error(sol_φ_approx,sol_φ)

## test NtD map without corners
D_D1,rhsmatrix = DF.obtain_ntd_matrices(Dmatrix,Smatrix,domain)
ev = DF.edge_indices(domain)    # edge variables
reduced_sol_φ = sol_φ[ev]
reduced_∂sol∂n_φ = ∂sol∂n_φ[ev]
rhs2 = rhsmatrix*reduced_∂sol∂n_φ
sol_φ_edge_approx = D_D1 \ rhs2
err_ntd_schur = DF.rel_error(sol_φ_edge_approx,reduced_sol_φ)

## try reduced indexing
bottom_idx = DF.reduced_boundary_indices(domain,:bottom)
left_idx = DF.reduced_boundary_indices(domain,:left)
right_idx = DF.reduced_boundary_indices(domain,:right) |> reverse  # REVERSED
top_idx = DF.reduced_boundary_indices(domain,:top)
@assert length(ev) == length(bottom_idx)+length(left_idx)+length(right_idx)+length(top_idx)

# check Bloch boundary conditions
utop = reduced_sol_φ[top_idx]
ubottom = reduced_sol_φ[bottom_idx]
∂utop = reduced_∂sol∂n_φ[top_idx]
∂ubottom = reduced_∂sol∂n_φ[bottom_idx]
v = reduced_sol_φ[left_idx]
∂v = reduced_∂sol∂n_φ[left_idx]
w = reduced_sol_φ[(right_idx)]   #[right_idx]
∂w = reduced_∂sol∂n_φ[(right_idx)]
norm_1 = norm(v*γ - w)
norm_2 = norm(∂v*γ + ∂w)

## test V matrices
# V matrices
V = D_D1 \ rhsmatrix #obtain_ntd_map(domain)
V11 = V[bottom_idx,bottom_idx]
V12 = V[bottom_idx,left_idx]
V13 = V[bottom_idx,right_idx]
V14 = V[bottom_idx,top_idx]
V21 = V[left_idx,bottom_idx]
V22 = V[left_idx,left_idx]
V23 = V[left_idx,right_idx]
V24 = V[left_idx,top_idx]
V31 = V[right_idx,bottom_idx]
V32 = V[right_idx,left_idx]
V33 = V[right_idx,right_idx]
V34 = V[right_idx,top_idx]
V41 = V[top_idx,bottom_idx]
V42 = V[top_idx,left_idx]
V43 = V[top_idx,right_idx]
V44 = V[top_idx,top_idx]

V_assemble = [V11 V12 V13 V14;
              V21 V22 V23 V24;
              V31 V32 V33 V34;
              V41 V42 V43 V44]
u_assemble = vcat(ubottom,v,w,utop)
∂u_assemble = vcat(∂ubottom,∂v,∂w,∂utop)

norm_v1 = norm(V*reduced_∂sol∂n_φ-reduced_sol_φ)
norm_v2 = norm(V_assemble*∂u_assemble-u_assemble)

## test N matrices
_V13 = -V13
_V23 = -V23
_V33 = -V33
_V43 = -V43
# C matrices
C1 = V12 + γ*_V13
C2 = V42 + γ*_V43
# D matrices
D0 = γ*V22 + γ^2*_V23 - V32 - γ*_V33
D1 = D0 \ (V31-γ*V21)
D2 = D0 \ (V34-γ*V24)
# N matrices
N11 = V11 + C1*D1
N12 = V14 + C1*D2
N21 = V41 + C2*D1
N22 = V44 + C2*D2

N_assemble = [N11 N12;
              N21 N22]
u_border = vcat(ubottom,utop)
∂u_border = vcat(∂ubottom,∂utop)

err_n = DF.rel_error(N_assemble*∂u_border,u_border)
#err_n_inv = DF.rel_error(N_assemble \ u_border,∂u_border)  # doesn't work, ill-conditioned

## test Q and Y matrices
# initial conditions
Q = -geo.Gbottom*(-1)     # REVERSE SIGN!!
Y = I(size(Q,1))
γ = exp(im*geo.α0*geo.L)
# generate data
top_qpoints = (DF.qpoint(domain,i) for i in DF.topboundary_indices(domain))
bottom_qpoints = (DF.qpoint(domain,i) for i in DF.bottomboundary_indices(domain))

u_top = [sol(DF.point(q)) for q in top_qpoints]
∂u_top = [∂sol∂n(DF.point(q),DF.normal(q)) for q in top_qpoints]
u_bottom = [sol(DF.point(q)) for q in bottom_qpoints]
∂u_bottom = [∂sol∂n(DF.point(q),DF.normal(q)) for q in bottom_qpoints]

# check 
err_q0 = DF.rel_error(Q*u_bottom,∂u_bottom)
err_q0_inv = DF.rel_error(Q\∂u_bottom,u_bottom) # doesn't work, ill-conditioned

# check 2
err_n1 = DF.rel_error(N11*∂u_bottom+N12*∂u_top,u_bottom)
err_n2 = DF.rel_error(N21*∂u_bottom+N22*∂u_top,u_top)

# iterate
Z = (I-N11*Q) \ N12
u_bottom_apprx = Z*∂u_top
err_z = DF.rel_error(u_bottom_apprx,u_bottom)

err_eq2 = DF.rel_error(N21*Q*u_bottom_apprx+N22*∂u_top,u_top)
eq2_matrix = N21*Q*Z+N22
err_eq2_inv = DF.rel_error(eq2_matrix\u_top,∂u_top)  # doesn't work, ill-conditioned

Qtop_inv = N22+N21*Q*Z
Qtop = inv(Qtop_inv)
Ytop = Y*Z*Qtop

# solve 
err_qinv = DF.rel_error(Qtop_inv*∂u_top,u_top)
err_q = DF.rel_error(Qtop*u_top,∂u_top) # doesn't work, ill-conditioned
err_y = DF.rel_error(Ytop*u_top,u_bottom) # works!

## solve for bottom
u_bottom_apprx_1 = Ytop*u_top
err_y_1 = DF.rel_error(u_bottom_apprx_1,u_bottom) # works! true equation

u_bottom_apprx_2 = Y*Z*∂u_top
err_y_2 = DF.rel_error(u_bottom_apprx_2,u_bottom) # works!

∂u_top_apprx = Qtop*u_top
err_∂u_top_apprx = DF.rel_error(∂u_top_apprx,∂u_top)
err_∂u_top_apprx_vec = abs.(∂u_top_apprx-∂u_top)/maximum(abs.(∂u_top))
u_bottom_apprx_3 = Y*Z*∂u_top_apprx
err_y_2 = DF.rel_error(u_bottom_apprx_3,u_bottom)


