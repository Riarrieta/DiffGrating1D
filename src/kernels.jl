
## Helmholtz

function double_layer_kernel(tx,ty,τ::QPoint,k)
    # actually, -2*D, where D is the double layer kernel
    τx,τy = point(τ)
    τx_p,τy_p = tvector(τ)
    tτnorm = sqrt((tx-τx)^2+(ty-τy)^2)
    L = im*k/2*(τy_p*(τx-tx)-τx_p*(τy-ty))*hankel1(k*tτnorm)/tτnorm
    return L
end
function double_layer_kernel(t::QPoint,τ::QPoint,k)
    tx,ty = point(t)
    return double_layer_kernel(tx,ty,τ,k)
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

function single_layer_kernel(tx,ty,τ::QPoint,k)
    # actually, 2*S, where S is the single layer kernel
    τx,τy = point(τ)
    τ_pnorm = tnorm(τ)
    tτnorm = sqrt((tx-τx)^2+(ty-τy)^2)
    M = im/2*hankel0(k*tτnorm)*τ_pnorm
    return M
end
function single_layer_kernel(t::QPoint,τ::QPoint,k)
    tx,ty = point(t)
    return single_layer_kernel(tx,ty,τ,k)
end
function single_layer_kernel_M1_M2(t::QPoint,τ::QPoint,k)
    # M1
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
    # correction of diagonal term, in case a change of variable is used
    ∂w = wprime(τ)
    is_change_of_variable_used = !isnan(∂w)
    M2_diag = M2_diag + (is_change_of_variable_used)*(2*log(∂w))
    M2 = (are_equal)*M2_diag + (!are_equal)*M2
    return M1,M2
end

## Corrections, for domains with corners
function double_layer_kernel_laplace_correction(tx,ty,τ::QPoint)
    τx,τy = point(τ)
    τx_p,τy_p = tvector(τ)
    τx_pp,τy_pp = ttvector(τ)
    τ_pnorm = tnorm(τ)
    tτnorm = sqrt((tx-τx)^2+(ty-τy)^2)
    H = 1/π*(τy_p*(tx-τx)-τx_p*(ty-τy))/tτnorm^2
    # diagonal term
    are_equal = iszero(tτnorm)
    Hdiag = 1/π*(τy_p*τx_pp-τx_p*τy_pp)/τ_pnorm^2
    H = (are_equal)*Hdiag + (!are_equal)*H
    return H
end
function double_layer_kernel_laplace_correction(t::QPoint,τ::QPoint)
    tx,ty = point(t)
    return double_layer_kernel_laplace_correction(tx,ty,τ)
end

## Laplace (for testing purposes)

function single_layer_kernel_laplace(tx,ty,τ::QPoint)
    τx,τy = point(τ)
    τ_pnorm = tnorm(τ)
    tτnorm = sqrt((tx-τx)^2+(ty-τy)^2)
    M = -1/(2π)*log(tτnorm)*τ_pnorm
    return M
end
function single_layer_kernel_laplace(t::QPoint,τ::QPoint)
    tx,ty = point(t)
    return single_layer_kernel_laplace(tx,ty,τ)
end

function double_layer_kernel_laplace(tx,ty,τ::QPoint)
    τx,τy = point(τ)
    τ_nx,τ_ny = normal(τ)
    τ_pnorm = tnorm(τ)
    tτnorm = sqrt((tx-τx)^2+(ty-τy)^2)
    L = -1/(2π)/(tτnorm)^2*rrdot(Point2D{Float64}(τx-tx,τy-ty),Point2D{Float64}(τ_nx,τ_ny))*τ_pnorm
    return L
end
function double_layer_kernel_laplace(t::QPoint,τ::QPoint)
    tx,ty = point(t)
    return double_layer_kernel_laplace(tx,ty,τ)
end

## Helmholtz (for testing purposes)

function single_layer_kernel_helmholtz(tx,ty,τ::QPoint,k)
    τx,τy = point(τ)
    τ_pnorm = tnorm(τ)
    tτnorm = sqrt((tx-τx)^2+(ty-τy)^2)
    M = im/4*hankel0(k*tτnorm)*τ_pnorm
    return M
end
function single_layer_kernel_helmholtz(t::QPoint,τ::QPoint,k)
    tx,ty = point(t)
    return single_layer_kernel_helmholtz(tx,ty,τ,k)
end

function double_layer_kernel_helmholtz(tx,ty,τ::QPoint,k)
    τx,τy = point(τ)
    τ_nx,τ_ny = normal(τ)
    τ_pnorm = tnorm(τ)
    tτnorm = sqrt((tx-τx)^2+(ty-τy)^2)
    L = -k*im/4*hankel1(k*tτnorm)/tτnorm*rdot(Point2D{Float64}(τx-tx,τy-ty),Point2D{Float64}(τ_nx,τ_ny))*τ_pnorm
    return L
end
function double_layer_kernel_helmholtz(t::QPoint,τ::QPoint,k)
    tx,ty = point(t)
    return double_layer_kernel_helmholtz(tx,ty,τ,k)
end