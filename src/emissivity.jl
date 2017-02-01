# This file contains the functions to compute the emissivities of bulks and multilayers

using Cubature # find package for integration

BulkOrMultiLayer = Union{Bulk,MultiLayer}


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

# Total emissivity

function bose_einstein(w,T)
    ħ*w/(exp(ħ*w/kb/T)-1.0)
end


function emissivity(struct :: BulkOrMultiLayer,T)
    e(u) = emissivity_w(struct, u*kb*T/ħ)*u^3/(exp(u)-1.0)
    e2(t) = e(t/(1.0-t))/(1.0-t)^2

    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = hquadrature(e2, 0.0, 1.0 ; reltol=1e-8, abstol=0, maxevals=0)

    return val*kb^4/ħ^3/c0^2/(2.0*pi)^2/sigma
end
