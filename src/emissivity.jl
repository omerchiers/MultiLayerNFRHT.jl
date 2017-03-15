# This file contains the functions to compute the emissivities of bulks and multilayers

using Cubature # find package for integration

BulkOrMultiLayer = Union{Layer,Bulk,MultiLayer}
LayerOrMultiLayer = Union{Layer,MultiLayer}

# Factors
function bose_einstein(w,T)
    ħ*w/(exp(ħ*w/kb/T)-1.0)
end

function planck(w,T)
    u = w*ħ/kb/T
    return  u^3/(exp(u)-1.0)*kb^4*T^4/ħ^3/c0^2/(2.0*pi)^2
end

# Wien frequency : return the frequency for which Planck distribution is maximum
function wien(T)
    return   2.8214393721220787*kb*T/ħ
end

function emissivity_kx_w(struct :: BulkOrMultiLayer, kx, w)
# Monocromatic emissivity

    (rte,t)=rt(struct, te(),kx,w)
    (rtm,t)=rt(struct, tm(),kx,w)

    integr1 = 1.0-abs(rte)^2
    integr2 = 1.0-abs(rtm)^2

    return integr1+integr2

end


function emissivity_w(struct :: BulkOrMultiLayer, w)
    e_kx(kx) = kx*emissivity_kx_w(struct ,kx, w)

    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = hquadrature(e_kx, 0.0, w/c0 ; reltol=1e-8, abstol=0, maxevals=0)

    return val/(w/c0)^2
end


function emissivity(struct :: BulkOrMultiLayer,T)
    e(u) = emissivity_w(struct, u*kb*T/ħ)*u^3/(exp(u)-1.0)
    e2(t) = e(t/(1.0-t))/(1.0-t)^2

    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = hquadrature(e2, 0.0, 1.0 ; reltol=1e-8, abstol=0, maxevals=0)

    return val*kb^4/ħ^3/c0^2/(2.0*pi)^2/sigma
end

function emissivity(struct :: BulkOrMultiLayer,T,wi,wf)
    e(u) = emissivity_w(struct, u*kb*T/ħ)*u^3/(exp(u)-1.0)
    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = hquadrature(e, wi*ħ/kb/T , wf*ħ/kb/T ; reltol=1e-8, abstol=0, maxevals=0)

    return val*kb^4/ħ^3/c0^2/(2.0*pi)^2/sigma
end
