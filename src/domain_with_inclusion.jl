struct DomainWithInclusion <: AbstractDomain
    ext_dom::DomainMultipleCorners
    int_dom::Domain
end
∂w(::QPoint,::Domain) = true  # change of variable is identity
∂w(q::QPoint,::AbstractDomainWithCorners) = wprime(q)

function double_layer_potential_matrix(ext_dom::AbstractDomain,domA::AbstractDomain,domB::AbstractDomain)
    # domA: integration domain
    # domB: observation domain
    # computed with trapezoidal rule
    NA = nunknowns(domA)
    na = domA.n
    ka = wavenumber(ext_dom)
    NB = nunknowns(domB)
    D = Array{ComplexF64}(undef, NB, NA)
    for j in edge_indices(domA) # skip corner nodes
        qj = qpoint(domA,j)
        ∂wj = ∂w(qj,domA)
        for i in 1:NB 
            qi = qpoint(domB,i)
            D[i,j] = -π/na*double_layer_kernel(qi,qj,ka)*∂wj
        end
    end
    # corner nodes
    for j in corner_indices(domA)
        for i in 1:NB  
            D[i,j] = zero(ComplexF64)
        end
    end 
    return D
end

function single_layer_potential_matrix(ext_dom::AbstractDomain,domA::AbstractDomain,domB::AbstractDomain)
    # domA: integration domain
    # domB: observation domain
    # computed with trapezoidal rule
    NA = nunknowns(domA)
    na = domA.n
    ka = wavenumber(ext_dom)
    NB = nunknowns(domB)
    S = Array{ComplexF64}(undef, NB, NA)
    for j in edge_indices(domA) # skip corner nodes
        qj = qpoint(domA,j)
        ∂wj = ∂w(qj,domA)
        for i in 1:NB 
            qi = qpoint(domB,i)
            S[i,j] = π/na*single_layer_kernel(qi,qj,ka)*∂wj
        end
    end
    # corner nodes
    for j in corner_indices(domA)
        for i in 1:NB  
            S[i,j] = zero(ComplexF64)
        end
    end 
    return S
end

function _assemble_T_G_matrices(d::DomainWithInclusion,boundary_label::Symbol,βvec,αlist,L)
    return _assemble_T_G_matrices(d.ext_dom,boundary_label::Symbol,βvec,αlist,L)
end

boundary_indices(d::DomainWithInclusion,label::Symbol) = boundary_indices(d.ext_dom,label)
reduced_boundary_indices(d::DomainWithInclusion,label::Symbol) = reduced_boundary_indices(d.ext_dom,label)
topboundary_indices(d::DomainWithInclusion) = topboundary_indices(d.ext_dom)
bottomboundary_indices(d::DomainWithInclusion) = bottomboundary_indices(d.ext_dom)
leftboundary_indices(d::DomainWithInclusion) = leftboundary_indices(d.ext_dom)
rightboundary_indices(d::DomainWithInclusion) = rightboundary_indices(d.ext_dom)
corner_indices(d::DomainWithInclusion) = corner_indices(d.ext_dom)
edge_indices(d::DomainWithInclusion) = edge_indices(d.ext_dom)
qpoint(d::DomainWithInclusion,i::Integer) = qpoint(d.ext_dom,i)

_check_ntd(d::DomainWithInclusion,domain_idx,α0,γ,ϵtol) = _check_ntd(d.ext_dom,domain_idx,α0,γ,ϵtol)

## Plot recipe
@recipe function plot_domain(d::DomainWithInclusion;tangent=true,normal=true)
    tangent := tangent
    normal := normal
    @series begin
        d.int_dom
    end
    return d.ext_dom
end