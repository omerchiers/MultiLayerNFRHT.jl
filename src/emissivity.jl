# This file contains the functions to compute the emissivities of bulks and multilayers

# Factors
"""
    bose_einstein(w,T)

Bose-Einstein factor using the definition:
```math
\\Theta(\\omega,T) = \\frac{\\hbar\\omega}{e^{\\frac{\\hbar\\omega}{k_B T}}-1}
```

"""
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
    planck(w,T)

Planck's Distribution using the definition:

```math
q^{\\text{BB}}(\\omega,T) =  \\Theta(\\omega,T) \\frac{k_0^2}{4\\pi^2}
```
"""
function planck(w,T)
    return  bose_einstein(w,T)*(w/c0)^2/(2.0*pi)^2
end

"""
    planck_fraction(w1,w2,T)

fraction of the energy between pulsations w1 and w2

```math
q^{\\text{BB}}(\\omega_1,\\omega_2,T) =  \\frac{1}{\\sigma T^4}\\int_{\\omega_1}^{\\omega_2}q^{\\text{BB}}(\\omega)d\\omega
```
"""
function planck_fraction(w1,w2,T)
  pl(w) = planck(w,T)
  (val, err) = quadgk(pl, w1, w2 ;rtol=1e-8)
  return val/sigma/T^4
end

"""
    wien(T)

Wien frequency: return the frequency for which Planck distribution is maximum.
T is given in Kelvin and the result is given in radHz
"""
function wien(T)
    return   2.8214393721220787*kb*T/ħ
end


"""
    invwien(w)

Inverse of Wien law : return the temperature (K)
for a given frequency in radHz
"""
function invwien(w)
    return   w*ħ/(2.8214393721220787*kb)
end

" Unit conversion from rad/s to wavelength in m"
unitconv(w) = 2.0*pi*c0/w

" Wien's wavelength obtained from wien's law in m"
lambda_wien(T) = 2.8977729e-3/T

"Monocromatic directional emissivity"
function emissivity_kx_w(structure, kx, w)

    (Rte,Tte)=power_rt(structure, te(), kx ,w)
    (Rtm,Ttm)=power_rt(structure, tm(), kx ,w)

    if imag(permittivity(substrate(structure),w)) == 0.0
        integr1 = 1.0 - Rte - Tte
        integr2 = 1.0 - Rtm - Ttm
    else
        integr1 = 1.0 - Rte
        integr2 = 1.0 - Rtm
    end

    return integr1 + integr2
end



"Total directional emissivity"
function emissivity_kx(structure,T, kx, wi,wf; rtol=1e-8)

    e(u) = emissivity_kx_w(structure, kx, u*kb*T/ħ)*u^3/(exp(u)-1.0)
    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = quadgk(e, wi*ħ/kb/T , wf*ħ/kb/T ; rtol=rtol)

    return val*kb^4/ħ^3/c0^2/(2.0*pi)^2/sigma
end

" Monocromatic hemispherical emissivity"
function emissivity_w(structure, w; rtol = 1e-8)
    e_kx(kx) = kx*emissivity_kx_w(structure ,kx, w)

    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = quadgk(e_kx, 0.0, w/c0 ; rtol=rtol)

    return val/(w/c0)^2
end

" total hemispherical emissivity"
function emissivity(structure,T;rtol = 1e-8)
    e(u) = emissivity_w(structure, u*kb*T/ħ;rtol = rtol)*u^3/(exp(u)-1.0)
    e2(t) = e(t/(1.0-t))/(1.0-t)^2

    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = quadgk(e2, 0.0, 1.0 ; rtol=rtol)

    return val*kb^4/ħ^3/c0^2/(2.0*pi)^2/sigma
end

" total hemispherical emissivity with integration bounds for frequency"
function emissivity(structure,T,wi,wf;rtol = 1e-8)
    e(u) = emissivity_w(structure, u*kb*T/ħ; rtol = rtol)*u^3/(exp(u)-1.0)
    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = quadgk(e, wi*ħ/kb/T , wf*ħ/kb/T ; rtol=rtol)

    return val*kb^4/ħ^3/c0^2/(2.0*pi)^2/sigma
end

function emissivity_fraction(structure,T,wi,wf;rtol = 1e-8)
    return emissivity_kx(structure,T,0.0,wi,wf;rtol = rtol)/planck_fraction(wi,wf,T)
end
