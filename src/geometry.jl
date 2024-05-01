
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
end

# check if a domain has a flat boundary (in the x-direction) of length L
function _check_flat_boundary(domain,boundary_label,L)
    idx = boundary_indices(domain,boundary_label)
    check1 = abs(qpoint(domain,idx[end]).x-qpoint(domain,idx[1]).x) ≤ L  # horizontal distance less than L
    check2 = all(iszero(qpoint(domain,i).ty) for i in idx)
    return check1 && check2
end

# check that bottom and top domain have the same boundary points
function _check_common_nodes(bottom_domain,top_domain)
    bottom_domain_top_nodes = (qpoint(bottom_domain,i) for i in boundary_indices(bottom_domain,:top))
    top_domain_bottom_nodes = (qpoint(top_domain,i) for i in boundary_indices(top_domain,:bottom))
    check = all(b.x ≈ t.x && b.y ≈ t.y for (b,t) in zip(bottom_domain_top_nodes,top_domain_bottom_nodes))
    return check
end
function _check_common_nodes(domains)
    check = true
    for i in 1:(length(domains)-1)
        bottom_domain = domains[i]
        top_domain = domains[i+1]
        _check = _check_common_nodes(bottom_domain,top_domain)
        _check || @info "" i
        check &= _check
    end
    return check
end

function Geometry(domains,kbottom,ktop,L,α0,Jmax)
    Jrange = -Jmax:Jmax
    αlist = [α0 + 2π*j/L for j in Jrange]
    βbottom = [sqrt(Complex(kbottom^2-aj^2)) for aj in αlist]
    βtop = [sqrt(Complex(ktop^2-aj^2)) for aj in αlist]
    # check top and bottom domains have flat boundaries 
    @assert _check_flat_boundary(domains[end],:top,L)
    @assert _check_flat_boundary(domains[1],:bottom,L)
    # check that domains have common boundaries
    @assert _check_common_nodes(domains)
    return Geometry(domains,kbottom,ktop,L,α0,Jmax,αlist,βbottom,βtop) 
end