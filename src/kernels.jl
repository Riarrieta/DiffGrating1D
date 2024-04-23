
function double_layer_kernel(t::QPoint,τ::QPoint,k)
    # actually, -2*D, where D is the double layer kernel
    tx,ty = point(t)
    τx,τy = point(τ)
    τx_p,τy_p = tvector(τ)
    tτnorm = distance(t,τ)
    L = im*k/2*(τy_p*(τx-tx)-τx_p*(τy-ty))*hankel1(k*tτnorm)/tτnorm
    return L
end
function double_layer_kernel_L1_L2(t::QPoint,τ::QPoint,k)
    # L1
    tx,ty = point(t)
    τx,τy = point(τ)
    τx_p,τy_p = tvector(τ)
    tτnorm = distance(t,τ)
    L1 = k/(2π)*(-τy_p*(τx-tx)+τx_p*(τy-ty))*besselj1(k*tτnorm)/tτnorm
    are_equal = iszero(tτnorm)
    L1 = (!are_equal)*L1 # diagonal term
    # L2
    t_param = param(t)
    τ_param = param(τ)
    L2 = double_layer_kernel(t,τ,k) - L1*log(4*sin((t_param-τ_param)/2)^2)
    # diagonal term
    tx_p,ty_p = tvector(t)
    tx_pp,ty_pp = ttvector(t)
    t_pnorm2 = tnorm2(t)
    L2diag = 1/(2π)*(tx_p*ty_pp-ty_p*tx_pp)/t_pnorm2
    L2 = (are_equal)*L2diag + (!are_equal)*L2
    return L1,L2
end

function single_layer_kernel(t::QPoint,τ::QPoint,k)
    # actually, -2*S, where S is the single layer kernel
    τ_pnorm = tnorm(τ)
    tτnorm = distance(t,τ)
    M = im/2*hankel0(k*tτnorm)*τ_pnorm
    return M
end
function single_layer_kernel_M1_M2(t::QPoint,τ::QPoint,k)
    # M1
    tx,ty = point(t)
    τx,τy = point(τ)
    τx_p,τy_p = tvector(τ)
    τ_pnorm = tnorm(τ)
    tτnorm = distance(t,τ)
    M1 = -1/(2π)*besselj0(k*tτnorm)*τ_pnorm
    # M2
    t_param = param(t)
    τ_param = param(τ)
    M2 = single_layer_kernel(t,τ,k) - M1*log(4*sin((t_param-τ_param)/2)^2)
    # diagonal term
    are_equal = iszero(tτnorm)
    C = Base.MathConstants.eulergamma
    t_pnorm = tnorm(t)
    M2_diag = (im/2-C/π-1/π*log(k/2*t_pnorm))*t_pnorm
    M2 = (are_equal)*M2_diag + (!are_equal)*M2
    return M1,M2
end

