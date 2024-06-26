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
    check_foo(q1,q2) = q1.x≈q2.x && q1.y≈q2.y && rdot(normal(q1),normal(q2))≈-1  # same points, opposite normals
    check = all(check_foo(b,t) for (b,t) in zip(bottom_domain_top_nodes,top_domain_bottom_nodes))
    #println(collect((q1.x,q2.x,q1.y,q2.y,rdot(normal(q1),normal(q2))) for (q1,q2) in zip(bottom_domain_top_nodes,top_domain_bottom_nodes)))
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
function _check_ntd(domain::AbstractDomain,domain_idx,α0,γ,ϵtol)
    k = wavenumber(domain)
    β = sqrt(Complex(k^2-α0^2))
    #x0 = Point2D(-5.0,-1.0)     # exterior eval point
    #sol(x) = hankel0(k*norm(x-x0))
    #sol_grad(x) = -k*hankel1(k*norm(x-x0))*(x-x0)/norm(x-x0)
    #w = k/sqrt(2)
    #sol(x) = cos(w*x[1])*cos(w*x[2])
    #sol_grad(x) = ForwardDiff.gradient(sol,x)
    sol(x) = exp(im*(α0*x[1] + β*x[2]))
    sol_grad(x) = Point2D(im*α0,im*β)*sol(x)
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
    D_D1,rhsmatrix = obtain_ntd_matrices(Dmatrix,Smatrix,domain)
    ev = edge_indices(domain)    # edge variables
    rhs2 = rhsmatrix*∂sol∂n_φ[ev]
    sol_φ_edge_approx = D_D1 \ rhs2
    err_ntd_schur = rel_error(sol_φ_edge_approx,sol_φ[ev])
    @info "NtD Schur" domain_idx err_ntd_schur
    check_ntd_schur = err_ntd_schur < ϵtol

    ## check N matrices
    err_n = _check_ntd_N_matrices(domain,sol_φ,∂sol∂n_φ,γ)
    @info "N matrices" domain_idx err_n
    check_n_matrices = err_n < ϵtol

    return check_ntd && check_ntd_schur && check_n_matrices
end
function _check_ntd(geo::Geometry,ϵtol)
    α0 = geo.α0
    γ = γfactor(geo)
    check = true
    for (i,domain) in enumerate(geo.domains)
        check &= _check_ntd(domain,i,α0,γ,ϵtol)
    end
    return check
end

function _check_ntd_N_matrices(domain::AbstractDomain,sol_φ,∂sol∂n_φ,γ)
    N11,N12,N21,N22 = obtain_reduced_ntd_map(domain,γ)
    u_bottom = sol_φ[bottomboundary_indices(domain)]
    u_top = sol_φ[topboundary_indices(domain)]
    ∂u_bottom = ∂sol∂n_φ[bottomboundary_indices(domain)]
    ∂u_top = ∂sol∂n_φ[topboundary_indices(domain)]
    u_bottom_approx = N11*∂u_bottom + N12*∂u_top
    u_top_approx = N21*∂u_bottom + N22*∂u_top
    err_n = max(rel_error(u_bottom_approx,u_bottom),rel_error(u_top_approx,u_top))
    return err_n
end

function check_geometry(geo::Geometry;ϵtol=1e-4)
    # check top and bottom domains have flat boundaries 
    check_flat_top = _check_flat_boundary(topdomain(geo),:top,geo.L)
    check_flat_bottom = _check_flat_boundary(bottomdomain(geo),:bottom,geo.L)
    # check that domains have common boundaries
    check_common = _check_common_nodes(geo.domains)
    # check matrices
    check_matrices = _check_T_G_matrices(geo,ϵtol)
    # check NtD 
    check_ntd = _check_ntd(geo,ϵtol)
    # check total
    check_total = check_flat_top && check_flat_bottom && check_common && check_matrices && check_ntd
    @info "Geometry tests" check_flat_top check_flat_bottom check_common check_matrices check_ntd check_total
    return check_total
end

function _is_homogeneous(geo::Geometry)
    k = geo.kbottom
    for dom in geo.domains
        if dom isa DomainWithInclusion
            (wavenumber(dom.int_dom)≠k || wavenumber(dom.ext_dom)≠k) && return false
        else
            wavenumber(dom)≠k && return false
        end
        
    end
    return true
end
function check_homogeneous_geometry(geo::Geometry;ϵtol=1e-4)
    @assert _is_homogeneous(geo) "Geometry is not homogeneous"
    # solve 'homogeneous' diffraction problem and compare with true solution
    u_reflected,r_coeff,u_transmitted,t_coeff = solve_diffraction_problem(geo);
    # true solution
    topdom = topdomain(geo)
    bottomdom = bottomdomain(geo)
    bottom_points = (point(qpoint(bottomdom,i)) for i in bottomboundary_indices(bottomdom))
    α0,β0 = αβfactors(geo)
    u_inc(x) = exp(im*(α0*x[1]-β0*x[2]))   # incident planewave
    u_transmitted_true = [u_inc(x) for x in bottom_points]   
    t_coeff_true = zeros(size(t_coeff)); t_coeff_true[geo.Jmax+1]=1;  # only zeroth order transmits
    # compare
    err_reflected = norm(u_reflected,Inf)   # should go to zero
    err_r_coeff = norm(r_coeff,Inf)   # should go to zero
    err_transmitted = rel_error(u_transmitted,u_transmitted_true)
    err_t_coeff = rel_error(t_coeff,t_coeff_true)
    check_hom = err_reflected<ϵtol && err_r_coeff<ϵtol && err_transmitted<ϵtol && err_t_coeff<ϵtol
    @info "Homogeneous Geometry tests" err_reflected err_r_coeff err_transmitted err_t_coeff check_hom
    return check_hom
end