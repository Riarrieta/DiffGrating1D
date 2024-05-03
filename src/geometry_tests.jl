# Useful scripts to test Geometry

# check if a domain has a flat boundary (in the x-direction) of length L
function _check_flat_boundary(domain,boundary_label::Symbol,L)
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

# check Tmatrix and Gmatrix
function _check_T_G_matrices(geo::Geometry,domain::AbstractDomain,boundary_label::Symbol,Tmatrix,Gmatrix,bvec,ϵtol)
    J0 = geo.Jmax
    Jrange = -J0:J0
    tvec = [rand(ComplexF64) for _ in Jrange]
    xvec = (qpoint(domain,i).x for i in boundary_indices(domain,boundary_label))
    f(x) = sum(tj*exp(im*αj*x) for (tj,αj) in zip(tvec,geo.αlist))
    fvec = f.(xvec)
    g(x) = sum(tj*im*bj*exp(im*αj*x) for (tj,bj,αj) in zip(tvec,bvec,geo.αlist))
    gvec = g.(xvec)

    tvec_approx = Tmatrix*fvec
    err_t = rel_error(tvec_approx,tvec)
    check_t = err_t < ϵtol
    @info "Tmatrix" boundary_label err_t

    gvec_approx = Gmatrix*fvec
    err_g = rel_error(gvec_approx,gvec)
    check_g = err_g < ϵtol
    @info "Gmatrix" boundary_label err_g
    return check_t && check_g
end
function _check_T_G_matrices(geo::Geometry,ϵtol)
    check_bottom = _check_T_G_matrices(geo,bottomdomain(geo),:bottom,geo.Tbottom,geo.Gbottom,geo.βbottom,ϵtol)
    check_top = _check_T_G_matrices(geo,topdomain(geo),:top,geo.Ttop,geo.Gtop,geo.βtop,ϵtol)
    return check_bottom && check_top
end

# tests NtD maps
function _check_ntd(domain::AbstractDomain,domain_idx,ϵtol)
    k = wavenumber(domain)
    x0 = Point2D(-5.0,-1.0)     # exterior eval point
    sol(x) = hankel0(k*norm(x-x0))
    sol_grad(x) = -k*hankel1(k*norm(x-x0))*(x-x0)/norm(x-x0)
    #sol(x) = cos(k*x[1])
    #sol_grad(x) = ForwardDiff.gradient(sol,x)
    ∂sol∂n(x,n) = rdot(sol_grad(x),n)
    sol_φ = [sol(point(q)) for q in domain.quad]
    ∂sol∂n_φ = [∂sol∂n(point(q),normal(q)) for q in domain.quad]

    ## matrices
    Dmatrix = double_layer_matrix_plus_identity(domain)
    Smatrix = single_layer_matrix(domain)
    #@info "" cond(Dmatrix)
    ## Ntd map
    sol_φ_approx = Dmatrix \ (Smatrix*∂sol∂n_φ)
    err_ntd = rel_error(sol_φ_approx,sol_φ)
    @info "NtD" domain_idx err_ntd
    check_ntd = err_ntd < ϵtol
    err_ntd_forward = rel_error(Dmatrix*sol_φ,Smatrix*∂sol∂n_φ)
    @info "NtD forward" domain_idx err_ntd_forward

    ## try Schur complement to eliminate corner variables
    cv = corner_indices(domain)    # corner variables
    ev = edge_indices(domain)    # edge variables
    D1 = Dmatrix[cv,cv]
    D2 = Dmatrix[cv,ev]
    D3 = Dmatrix[ev,cv]
    D4 = Dmatrix[ev,ev]
    S2 = Smatrix[cv,ev]
    S4 = Smatrix[ev,ev]

    D_D1 = (D4-D3*(D1\D2))  # Schur complement of D1
    rhsmatrix = (S4-D3*(D1\S2))
    rhs2 = rhsmatrix*∂sol∂n_φ[ev]
    sol_φ_edge_approx = D_D1 \ rhs2
    err_ntd_schur = rel_error(sol_φ_edge_approx,sol_φ[ev])
    @info "NtD Schur" domain_idx err_ntd_schur
    check_ntd_schur = err_ntd_schur < ϵtol
    return check_ntd && check_ntd_schur
end
function _check_ntd(geo::Geometry,ϵtol)
    check = true
    for (i,domain) in enumerate(geo.domains)
        check &= _check_ntd(domain,i,ϵtol)
    end
    return check
end

function check_geometry(geo::Geometry;ϵtol=1e-4)
    # check top and bottom domains have flat boundaries 
    check_boundary = _check_flat_boundary(topdomain(geo),:top,geo.L)
    check_flat = _check_flat_boundary(bottomdomain(geo),:bottom,geo.L)
    # check that domains have common boundaries
    check_common = _check_common_nodes(geo.domains)
    # check matrices
    check_matrices = _check_T_G_matrices(geo,ϵtol)
    # check NtD 
    check_ntd = _check_ntd(geo,ϵtol)
    # check total
    check_total = check_boundary && check_flat && check_common && check_matrices && check_ntd
    @info "Geometry tests" check_boundary check_flat check_common check_matrices check_ntd check_total
    return check_total
end