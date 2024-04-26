
## Potentials
function double_layer_potential(d::Domain,ϕ)
    # computed with trapezoidal rule
    function pot(x)
        result = zero(ComplexF64)
        for (ϕi,qi) in zip(ϕ,d.quad)
            result += double_layer_kernel(x...,qi,d.k)*ϕi
        end
        result = (-0.5)*π/d.n*result # (-0.5) factor to correct the kernel
        return result
    end
    return pot
end
function double_layer_potential(d::DomainWith1Corner,ϕ)
    # computed with trapezoidal rule
    function pot(x)
        result = zero(ComplexF64)
        for i in 2:nunknowns(d)  # skip corner, since w'(corner) = 0
            ϕi = ϕ[i]
            qi = d.quad[i]
            ∂wi = d.wprime[i]
            result += double_layer_kernel(x...,qi,d.k)*ϕi*∂wi
        end
        result = (-0.5)*π/d.n*result # (-0.5) factor to correct the kernel
        return result
    end
    return pot
end

function single_layer_potential(d::Domain,ϕ)
    # computed with trapezoidal rule
    function pot(x)
        result = zero(ComplexF64)
        for (ϕi,qi) in zip(ϕ,d.quad)
            result += single_layer_kernel(x...,qi,d.k)*ϕi
        end
        result = (0.5)*π/d.n*result # (0.5) factor to correct the kernel
        return result
    end
    return pot
end
function single_layer_potential(d::DomainWith1Corner,ϕ)
    # computed with trapezoidal rule
    function pot(x)
        result = zero(ComplexF64)
        for i in 2:nunknowns(d)  # skip corner, since w'(corner) = 0
            ϕi = ϕ[i]
            qi = d.quad[i]
            ∂wi = d.wprime[i]
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
function mk_weight(qi::QPoint,qj::QPoint,d::Domain)
    return mk_weight(param(qi),param(qj),d.n)
end

function double_layer_matrix_plus_identity(d::Domain)
    N = nunknowns(d)
    n = d.n
    k = wavenumber(d)
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
function double_layer_matrix_plus_identity(d::DomainWith1Corner)
    N = nunknowns(d)
    n = d.n
    k = wavenumber(d)
    D = Array{ComplexF64}(undef, N, N)
    for j in 2:N # skip corner node j=1
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
    # corner node j=1
    for i in 1:N  
        qi = qpoint(d,i)
        D[i,1] = -one(ComplexF64)
        for l in 2:N # skip corner node
            ql = qpoint(d,l)
            ∂wl = wprime(ql)
            H = double_layer_kernel_laplace_correction(qi,ql)
            D[i,1] += -π/n*H*∂wl
        end
        # add identity on diagonal
        D[i,1] = D[i,1] + (i==1)
    end
    return D
end

function single_layer_matrix(d::Domain)
    N = nunknowns(d)
    n = d.n
    k = wavenumber(d)
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
function single_layer_matrix(d::DomainWith1Corner)
    N = nunknowns(d)
    n = d.n
    k = wavenumber(d)
    D = Array{ComplexF64}(undef, N, N)
    for j in 2:N # skip corner node j=1
        qj = qpoint(d,j)
        ∂wj = wprime(qj)
        for i in 1:N  
            qi = qpoint(d,i)
            M1,M2 = single_layer_kernel_M1_M2(qi,qj,k)
            D[i,j] = (mk_weight(qi,qj,d)*M1 + π/n*M2)*∂wj
        end
    end
    # corner node j=1
    for i in 1:N  
        D[i,1] = zero(ComplexF64)
    end
    return D
end