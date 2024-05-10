
## Potentials
function double_layer_potential(d::Domain,ϕ;k=wavenumber(d))
    # computed with trapezoidal rule
    function pot(x)
        result = zero(ComplexF64)
        for (ϕi,qi) in zip(ϕ,d.quad)
            result += double_layer_kernel(x...,qi,k)*ϕi
        end
        result = (-0.5)*π/d.n*result # (-0.5) factor to correct the kernel
        return result
    end
    return pot
end
function double_layer_potential(d::AbstractDomainWithCorners,ϕ)
    # computed with trapezoidal rule
    function pot(x)
        result = zero(ComplexF64)
        for i in edge_indices(d)  # skip corners, since w'(corner) = 0
            ϕi = ϕ[i]
            qi = d.quad[i]
            ∂wi = wprime(qi)
            result += double_layer_kernel(x...,qi,d.k)*ϕi*∂wi
        end
        result = (-0.5)*π/d.n*result # (-0.5) factor to correct the kernel
        return result
    end
    return pot
end

function single_layer_potential(d::Domain,ϕ;k=wavenumber(d))
    # computed with trapezoidal rule
    function pot(x)
        result = zero(ComplexF64)
        for (ϕi,qi) in zip(ϕ,d.quad)
            result += single_layer_kernel(x...,qi,k)*ϕi
        end
        result = (0.5)*π/d.n*result # (0.5) factor to correct the kernel
        return result
    end
    return pot
end
function single_layer_potential(d::AbstractDomainWithCorners,ϕ)
    # computed with trapezoidal rule
    function pot(x)
        result = zero(ComplexF64)
        for i in edge_indices(d)  # skip corner, since w'(corner) = 0
            ϕi = ϕ[i]
            qi = d.quad[i]
            ∂wi = wprime(qi)
            result += single_layer_kernel(x...,qi,d.k)*ϕi*∂wi
        end
        result = (0.5)*π/d.n*result # (0.5) factor to correct the kernel
        return result
    end
    return pot
end

## Operators

function mk_weight(ti,tj,n)
    # Martensen and Kussmaul weights R_j(ti)^(n)
    d = ti-tj
    return -2π/n*sum(1/m*cos(m*d) for m in 1:n-1) - π/n^2*cos(n*d)
end
function mk_weight(qi::QPoint,qj::QPoint,d::AbstractDomain)
    return mk_weight(param(qi),param(qj),d.n)
end

function double_layer_matrix_plus_identity(d::Domain;k=wavenumber(d))
    N = nunknowns(d)
    n = d.n
    D = Array{ComplexF64}(undef, N, N)
    for j in 1:N
        qj = qpoint(d,j)
        for i in 1:N
            qi = qpoint(d,i)
            L1,L2 = double_layer_kernel_L1_L2(qi,qj,k)
            D[i,j] = -(mk_weight(qi,qj,d)*L1 + π/n*L2)
            # add identity on diagonal
            D[i,j] = D[i,j] + (i==j)
        end
    end
    return D
end

function corner_modification_function(qi::QPoint,corner_index_l::Int64,d::AbstractDomainWithCorners)
    dist_ij(j) = distance(qi,qpoint(d,j))
    dist_lj(j) = distance(qpoint(d,corner_index_l),qpoint(d,j))
    foo(j) = ifelse(j==corner_index_l,one(ComplexF64),dist_ij(j)/dist_lj(j))
    return prod(foo(j) for j in corner_indices(d))
end
function corner_modification_function(::QPoint,::Int64,::DomainWith1Corner)
    return one(ComplexF64)
end
function double_layer_matrix_plus_identity(d::AbstractDomainWithCorners)
    N = nunknowns(d)
    n = d.n
    k = wavenumber(d)
    D = Array{ComplexF64}(undef, N, N)
    for j in edge_indices(d) # skip corner nodes
        qj = qpoint(d,j)
        ∂wj = wprime(qj)
        for i in 1:N  
            qi = qpoint(d,i)
            L1,L2 = double_layer_kernel_L1_L2(qi,qj,k)
            D[i,j] = -(mk_weight(qi,qj,d)*L1 + π/n*L2)*∂wj
            # add identity on diagonal
            D[i,j] = D[i,j] + (i==j)
        end
    end
    # corner nodes
    for j in corner_indices(d)
        for i in 1:N  
            qi = qpoint(d,i)
            D[i,j] = -one(ComplexF64)
            for l in edge_indices(d) # skip corner nodes
                ql = qpoint(d,l)
                ∂wl = wprime(ql)
                H = double_layer_kernel_laplace_correction(qi,ql)
                D[i,j] += -π/n*H*∂wl
            end
            # include correction function
            D[i,j] = D[i,j]*corner_modification_function(qi,j,d)
            # add identity on diagonal
            D[i,j] = D[i,j] + (i==j)
        end
    end 
    return D
end

function single_layer_matrix(d::Domain;k=wavenumber(d))
    N = nunknowns(d)
    n = d.n
    D = Array{ComplexF64}(undef, N, N)
    for j in 1:N
        qj = qpoint(d,j)
        for i in 1:N
            qi = qpoint(d,i)
            M1,M2 = single_layer_kernel_M1_M2(qi,qj,k)
            D[i,j] = mk_weight(qi,qj,d)*M1 + π/n*M2
        end
    end
    return D
end
function single_layer_matrix(d::AbstractDomainWithCorners)
    N = nunknowns(d)
    n = d.n
    k = wavenumber(d)
    D = Array{ComplexF64}(undef, N, N)
    for j in edge_indices(d) # skip corner nodes
        qj = qpoint(d,j)
        ∂wj = wprime(qj)
        for i in 1:N  
            qi = qpoint(d,i)
            M1,M2 = single_layer_kernel_M1_M2(qi,qj,k)
            D[i,j] = (mk_weight(qi,qj,d)*M1 + π/n*M2)*∂wj
        end
    end
    # corner nodes
    for j in corner_indices(d)
        for i in 1:N  
            D[i,j] = zero(ComplexF64)
        end
    end
    return D
end