
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
    k = wavenumber(d)
    D = Array{ComplexF64}(undef, N, N)
    for j in 1:N
        qj = qpoint(d,j)
        for i in 1:N
            qi = qpoint(d,i)
            L1,L2 = double_layer_kernel_L1_L2(qi,qj,k)
            D[i,j] = mk_weight(qi,qj,d)*L1 + π/n*L2
            # add identity on diagonal
            D[i,j] = D[i,j] + (i==j)
        end
    end
    return D
end

function single_layer_matrix(d::Domain)
    N = nunknowns(d)
    k = wavenumber(d)
    D = Array{ComplexF64}(undef, N, N)
    for j in 1:N
        qj = qpoint(d,j)
        for i in 1:N
            qi = qpoint(d,i)
            M1,M2 = double_layer_kernel_M1_M2(qi,qj,k)
            D[i,j] = mk_weight(qi,qj,d)*M1 + π/n*M2
        end
    end
    return D
end