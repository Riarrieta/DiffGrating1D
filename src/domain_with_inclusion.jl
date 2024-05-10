struct DomainWithInclusion
    ext_dom::DomainMultipleCorners
    int_dom::Domain
end
∂w(::QPoint,::Domain) = true  # change of variable is identity
∂w(q::QPoint,::AbstractDomainWithCorners) = wprime(q)

function double_layer_potential_matrix(domA::AbstractDomain,domB::AbstractDomain)
    # domA: integration domain
    # domB: observation domain
    # computed with trapezoidal rule
    NA = nunknowns(domA)
    na = domA.n
    ka = wavenumber(domA)
    NB = nunknowns(domB)
    D = Array{ComplexF64}(undef, NB, NA)
    for j in edge_indices(domA) # skip corner nodes
        qj = qpoint(domA,j)
        ∂wj = ∂w(qj,domA)
        for i in 1:NB 
            qi = qpoint(domB,i)
            D[i,j] = π/na*double_layer_kernel(qi,qj,ka)*∂wj
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

function single_layer_potential_matrix(domA::AbstractDomain,domB::AbstractDomain)
    # domA: integration domain
    # domB: observation domain
    # computed with trapezoidal rule
    NA = nunknowns(domA)
    na = domA.n
    ka = wavenumber(domA)
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
