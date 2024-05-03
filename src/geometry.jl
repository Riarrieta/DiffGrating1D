
struct Geometry
    domains::Vector{<:AbstractDomain}  # domains, from bottom to top
    kbottom::Float64              # wavenumber bottom domain
    ktop::Float64                 # wavenumber top domain
    L::Float64       # period
    α0::Float64      # Bloch wavenumber
    Jmax::Int64      # maximum Fourier mode
    αlist::Vector{Float64}      # αj wavevectors
    βbottom::Vector{ComplexF64}   # βj wavevectors, bottom
    βtop::Vector{ComplexF64}      # βj wavevectors, top
    # matrices
    Tbottom::Matrix{ComplexF64}   # from grid to Fourier components, bottom
    Gbottom::Matrix{ComplexF64}   # iB^(2) operator, grid to grid, bottom
    Ttop::Matrix{ComplexF64}      # from grid to Fourier components, top
    Gtop::Matrix{ComplexF64}      # iB^(1) operator, grid to grid, top
end
bottomdomain(g::Geometry) = g.domains[1]
topdomain(g::Geometry) = g.domains[end]

function _assemble_T_G_matrices(domain::AbstractDomain,boundary_label::Symbol,βvec,αlist,L)
    qpoints = (qpoint(domain,i) for i in boundary_indices(domain,boundary_label))
    xlist = (q.x for q in qpoints)
    ∂wlist = (wprime(q) for q in qpoints)
    tnorm_vec = (tnorm(q) for q in qpoints)  # tangent vectors norm
    Δs = 2π/nunknowns(domain)    # spacing in (global) parameter space
    # from grid to Fourier
    Tmatrix = [Δs/L*exp(-im*αi*xj)*tnj*∂wj for αi in αlist, (xj,∂wj,tnj) in zip(xlist,∂wlist,tnorm_vec)]  
    # from Fourier to grid
    Fmatrix = [exp(im*αj*xi) for xi in xlist, αj in αlist]
    # iB operator, grid to grid
    Gmatrix = Fmatrix*diagm(im*βvec)*Tmatrix
    return Tmatrix,Gmatrix
end

function Geometry(domains,kbottom,ktop,L,α0,Jmax)
    Jrange = -Jmax:Jmax
    αlist = [α0 + 2π*j/L for j in Jrange]
    βbottom = [sqrt(Complex(kbottom^2-aj^2)) for aj in αlist]
    βtop = [sqrt(Complex(ktop^2-aj^2)) for aj in αlist]
    # assemble top and bottom matrices
    Tbottom,Gbottom = _assemble_T_G_matrices(domains[1],:bottom,βbottom,αlist,L)
    Ttop,Gtop = _assemble_T_G_matrices(domains[end],:top,βtop,αlist,L)
    return Geometry(domains,kbottom,ktop,L,α0,Jmax,αlist,βbottom,βtop,Tbottom,Gbottom,Ttop,Gtop) 
end