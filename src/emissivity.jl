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
q_{\\omega}^{\\text{BB}} =  \\Theta(\\omega,T) \\frac{k_0^2}{4π^2}
```
"""
function planck(w,T)
    return  bose_einstein(w,T)*(w/c0)^2/(2.0*pi)^2
end

" Computes the fraction of the blackbody spectrum "
function planck_fraction(w1,w2,T)
  pl(w) = planck(w,T)
  (val, err) = hquadrature(pl, w1, w2 ;reltol=1e-8, abstol=0, maxevals=0)
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


function emissivity_kx_w(structure :: BulkOrMultiLayer, kx, w)

    (rte,t)=rt(structure, te(),kx,w)
    (rtm,t)=rt(structure, tm(),kx,w)

    integr1 = 1.0-abs(rte)^2
    integr2 = 1.0-abs(rtm)^2

    return integr1+integr2

end

" Monocromatic emissivity"
function emissivity_w(structure :: BulkOrMultiLayer, w)
    e_kx(kx) = kx*emissivity_kx_w(structure ,kx, w)

    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = hquadrature(e_kx, 0.0, w/c0 ; reltol=1e-8, abstol=0, maxevals=0)

    return val/(w/c0)^2
end


function emissivity(structure :: BulkOrMultiLayer,T)
    e(u) = emissivity_w(structure, u*kb*T/ħ)*u^3/(exp(u)-1.0)
    e2(t) = e(t/(1.0-t))/(1.0-t)^2

    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = hquadrature(e2, 0.0, 1.0 ; reltol=1e-8, abstol=0, maxevals=0)

    return val*kb^4/ħ^3/c0^2/(2.0*pi)^2/sigma
end

function emissivity(structure :: BulkOrMultiLayer,T,wi,wf)
    e(u) = emissivity_w(structure, u*kb*T/ħ)*u^3/(exp(u)-1.0)
    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = hquadrature(e, wi*ħ/kb/T , wf*ħ/kb/T ; reltol=1e-8, abstol=0, maxevals=0)

    return val*kb^4/ħ^3/c0^2/(2.0*pi)^2/sigma
end
