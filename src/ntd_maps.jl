
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
    for dom in geo.domains  # from bottom to top
        # multiply by -1, since the normals of adyacent domains point in opposite directions
        rmul!(Q,-one(ComplexF64))  
        # obtain Q and Y from Domain
        Q,Y = obtain_Q_Y_matrices(dom,γ,Q,Y)
    end
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
    return u_reflected,r_coeff,u_transmitted,t_coeff
end