
function obtain_ntd_matrices(Dmatrix,Smatrix,domain::AbstractDomain)
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
    return D_D1,rhsmatrix
end
function obtain_ntd_matrices(domain::AbstractDomain)
    ## matrices
    Dmatrix = double_layer_matrix_plus_identity(domain)
    Smatrix = single_layer_matrix(domain)
    return obtain_ntd_matrices(Dmatrix,Smatrix,domain)
end
function obtain_ntd_map(domain::AbstractDomain)
    D_D1,rhsmatrix = obtain_ntd_matrices(domain::AbstractDomain)
    return D_D1 \ rhsmatrix
end

function obtain_reduced_ntd_map(domain::AbstractDomain,γ)
    # right indices MUST be reversed, so that
    # f[right_idx] = γ*f[left_idx]  (Bloch boundary condition)
    # where f is either u or ∂u/∂n
    bottom_idx = reduced_boundary_indices(domain,:bottom)
    left_idx = reduced_boundary_indices(domain,:left)
    right_idx = reduced_boundary_indices(domain,:right) |> reverse
    top_idx = reduced_boundary_indices(domain,:top)
    # V matrices
    # the Vx3 matrices MUST have a minus sign (the ones affecting ∂u[right]),
    # because in my convention the normal vector always points outward,
    # but in the paper convention it always points to the right
    V = obtain_ntd_map(domain)
    V11 = V[bottom_idx,bottom_idx]
    V12 = V[bottom_idx,left_idx]
    V13 = -V[bottom_idx,right_idx]
    V14 = V[bottom_idx,top_idx]
    V21 = V[left_idx,bottom_idx]
    V22 = V[left_idx,left_idx]
    V23 = -V[left_idx,right_idx]
    V24 = V[left_idx,top_idx]
    V31 = V[right_idx,bottom_idx]
    V32 = V[right_idx,left_idx]
    V33 = -V[right_idx,right_idx]
    V34 = V[right_idx,top_idx]
    V41 = V[top_idx,bottom_idx]
    V42 = V[top_idx,left_idx]
    V43 = -V[top_idx,right_idx]
    V44 = V[top_idx,top_idx]
    # C matrices
    C1 = V12 + γ*V13
    C2 = V42 + γ*V43
    # D matrices
    D0 = γ*V22 + γ^2*V23 - V32 - γ*V33
    D1 = D0 \ (V31-γ*V21)
    D2 = D0 \ (V34-γ*V24)
    # N matrices
    N11 = V11 + C1*D1
    N12 = V14 + C1*D2
    N21 = V41 + C2*D1
    N22 = V44 + C2*D2
    #@info "" size(C1) size(C2) size(D0) size(D1) size(D2) size(N11) size(N12) size(N21) size(N22)
    return N11,N12,N21,N22
end

function obtain_Q_Y_matrices(domain::AbstractDomain,γ,Qprev,Yprev)
    N11,N12,N21,N22 = obtain_reduced_ntd_map(domain,γ)
    #@info "" size(N11) size(Qprev) size(N12)
    Z = (I-N11*Qprev) \ N12
    Q = inv(N22+N21*Qprev*Z)
    Y = Yprev*Z*Q
    return Q,Y
end

function solve_diffraction_problem(geo::Geometry)
    # initial conditions
    γ = γfactor(geo)
    α0,β0 = αβfactors(geo)
    Q = -geo.Gbottom  # -iB^(2)
    Y = I(size(Q,1))
    # iterate
    for (i,dom) in enumerate(geo.domains) # from bottom to top
        @info "Solving Domain $i"
        # multiply by -1, since the normals of adyacent domains point in opposite directions
        rmul!(Q,-one(ComplexF64))  
        # obtain Q and Y from Domain
        Q,Y = obtain_Q_Y_matrices(dom,γ,Q,Y)
    end
    @info "Solving Top Domain"
    # solve top boundary
    topdom = topdomain(geo)
    top_points = (point(qpoint(topdom,i)) for i in topboundary_indices(topdom))
    u_inc(x) = exp(im*(α0*x[1]-β0*x[2]))   # incident planewave
    u_incident = [u_inc(x) for x in top_points]   
    β0 = geo.βtop[geo.Jmax+1]      # β of incident planewave
    rhs = [-2*im*β0*u for u in u_incident]
    Lmatrix = Q-geo.Gtop  # Q-iB^(1)
    # solve and get fields
    utop = Lmatrix\rhs
    u_reflected = utop-u_incident
    u_transmitted = Y*utop
    # get reflection and transmission coeff
    r_coeff = geo.Ttop*u_reflected
    t_coeff = geo.Tbottom*u_transmitted
    @info "Done"
    return u_reflected,r_coeff,u_transmitted,t_coeff
end

## Domains with inclusion
function obtain_pure_ntd_matrices(d::DomainWithInclusion)
    ext_dom = d.ext_dom
    int_dom = d.int_dom
    # exterior domain
    DmatrixA = double_layer_matrix_plus_identity(ext_dom)
    SmatrixA = single_layer_matrix(ext_dom)
    # interior domain
    DmatrixB = double_layer_matrix_plus_identity(int_dom)
    SmatrixB = single_layer_matrix(int_dom)
    N2 = DmatrixB \ SmatrixB  # NtD of interior domain
    # interior domain, but k of exterior domain
    DmatrixB_kA = double_layer_matrix_plus_identity(int_dom;k=ext_dom.k)
    SmatrixB_kA = single_layer_matrix(int_dom;k=ext_dom.k)
    # layer potential matrices
    DAB = double_layer_potential_matrix(ext_dom,ext_dom,int_dom)
    SAB = single_layer_potential_matrix(ext_dom,ext_dom,int_dom)
    DBA = double_layer_potential_matrix(ext_dom,int_dom,ext_dom)
    SBA = single_layer_potential_matrix(ext_dom,int_dom,ext_dom)
    # auxiliary matrices
    lhs_matrix_1 = 2*N2+SmatrixB_kA-DmatrixB_kA*N2
    Z1 = lhs_matrix_1 \ SAB
    Z2 = lhs_matrix_1 \ DAB
    NZ1 = N2*Z1
    NZ2 = N2*Z2
    lhs_matrix_2 = DmatrixA - SBA*Z2 + DBA*NZ2
    rhs_matrix_2 = SmatrixA - SBA*Z1 + DBA*NZ1
    return lhs_matrix_2, rhs_matrix_2
end
function obtain_pure_ntd_map(d::DomainWithInclusion)
    lhs_matrix,rhs_matrix = obtain_pure_ntd_matrices(d)
    N1 = lhs_matrix \ rhs_matrix   # NtD of exterior domain
    return N1
end

function obtain_ntd_matrices(domain::DomainWithInclusion)
    ## matrices
    lhs_matrix,rhs_matrix = obtain_pure_ntd_matrices(domain)
    return obtain_ntd_matrices(lhs_matrix,rhs_matrix,domain)
end

# Same as solve_diffraction_problem, but one single Domain
# gets repeated vertically 'nlevels' times.
function solve_diffraction_for_multiple_equal_domains(geo::Geometry;nlevels)
    @assert length(geo.domains) == 1
    @assert geo.domains[1] isa DomainWithInclusion
    domain = geo.domains[1]
    # initial conditions
    γ = γfactor(geo)
    α0,β0 = αβfactors(geo)
    Q = -geo.Gbottom  # -iB^(2)
    Y = I(size(Q,1))
    # domain N matrices
    N11,N12,N21,N22 = obtain_reduced_ntd_map(domain,γ)
    # iterate
    for (i) in 1:nlevels # from bottom to top
        @info "Solving Domain $i"
        # multiply by -1, since the normals of adyacent domains point in opposite directions
        rmul!(Q,-one(ComplexF64))  
        # obtain Q and Y from Domain
        Z = (I-N11*Q) \ N12
        Q = inv(N22+N21*Q*Z)
        Y = Y*Z*Q
    end
    @info "Solving Top Domain"
    # solve top boundary
    topdom = topdomain(geo)
    top_points = (point(qpoint(topdom,i)) for i in topboundary_indices(topdom))
    L = geo.L
    Δz = (nlevels-1)*L
    u_inc(x) = exp(im*(α0*x[1]-β0*(x[2]+Δz)))   # incident planewave
    u_incident = [u_inc(x) for x in top_points]   
    β0 = geo.βtop[geo.Jmax+1]      # β of incident planewave
    rhs = [-2*im*β0*u for u in u_incident]
    Lmatrix = Q-geo.Gtop  # Q-iB^(1)
    # solve and get fields
    utop = Lmatrix\rhs
    u_reflected = utop-u_incident
    u_transmitted = Y*utop
    # get reflection and transmission coeff
    r_coeff = geo.Ttop*u_reflected
    t_coeff = geo.Tbottom*u_transmitted
    @info "Done"
    return u_reflected,r_coeff,u_transmitted,t_coeff
end