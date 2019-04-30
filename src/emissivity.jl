# This file contains the functions to compute the emissivities of bulks and multilayers


BulkOrMultiLayer = Union{Layer,Bulk,MultiLayer}
LayerOrMultiLayer = Union{Layer,MultiLayer}

# Factors
function bose_einstein(w,T)
    u = w*ħ/kb/T
    if T==0
        return 0.0
    end
    if w==0.0
        return kb*T
    end

    return u*kb*T/(exp(u)-1.0)
end

"""
Planck's Distribution using the definition :
```math
q_{\\omega}^{\\text{BB}} =  \\Theta(\\omega,T) \\frac{k_0^2}{4\\pi^2}
```
"""
function planck(w,T)
    return  bose_einstein(w,T)*(w/c0)^2/(2.0*pi)^2
end

" Computes the fraction of the blackbody spectrum "
function planck_fraction(w1,w2,T)
  pl(w) = planck(w,T)
  (val, err) = quadgk(pl, w1, w2 ;rtol=1e-8)
  return val/sigma/T^4
end

" Wien frequency : return the frequency for which Planck distribution is maximum"
function wien(T)
    return   2.8214393721220787*kb*T/ħ
end

" Inverse of Wien law : return the temperature for a given frequency"
function invwien(w)
    return   w*ħ/(2.8214393721220787*kb)
end

" Unit conversion from rad/s to wavelength in m"
unitconv(w) = 2.0*pi*c0/w

" Wien's wavelength obtained from wien's law in m"
lambda_wien(T) = 2.8977729e-3/T

"Monocromatic directional emissivity"
function emissivity_kx_w(structure :: BulkOrMultiLayer, kx, w)

    (rte,t)=rt(structure, te(),kx,w)
    (rtm,t)=rt(structure, tm(),kx,w)

    integr1 = 1.0-abs(rte)^2
    integr2 = 1.0-abs(rtm)^2

    return integr1+integr2
end


"Total directional emissivity"
function emissivity_kx(structure :: BulkOrMultiLayer,T, kx, wi,wf)

    e(u) = emissivity_kx_w(structure, kx, u*kb*T/ħ)*u^3/(exp(u)-1.0)
    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = quadgk(e, wi*ħ/kb/T , wf*ħ/kb/T ; rtol=1e-8)

    return val*kb^4/ħ^3/c0^2/(2.0*pi)^2/sigma
end

" Monocromatic hemispherical emissivity"
function emissivity_w(structure :: BulkOrMultiLayer, w)
    e_kx(kx) = kx*emissivity_kx_w(structure ,kx, w)

    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = quadgk(e_kx, 0.0, w/c0 ; rtol=1e-8)

    return val/(w/c0)^2
end

" total hemispherical emissivity"
function emissivity(structure :: BulkOrMultiLayer,T)
    e(u) = emissivity_w(structure, u*kb*T/ħ)*u^3/(exp(u)-1.0)
    e2(t) = e(t/(1.0-t))/(1.0-t)^2

    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = quadgk(e2, 0.0, 1.0 ; rtol=1e-8)

    return val*kb^4/ħ^3/c0^2/(2.0*pi)^2/sigma
end

" total hemispherical emissivity with integration bounds for frequency"
function emissivity(structure :: BulkOrMultiLayer,T,wi,wf)
    e(u) = emissivity_w(structure, u*kb*T/ħ)*u^3/(exp(u)-1.0)
    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = quadgk(e, wi*ħ/kb/T , wf*ħ/kb/T ; rtol=1e-8)

    return val*kb^4/ħ^3/c0^2/(2.0*pi)^2/sigma
end

function emissivity_fraction(structure :: BulkOrMultiLayer,T,wi,wf)
    return emissivity_kx(structure,T,0.0,wi,wf)/planck_fraction(wi,wf,T)
end
