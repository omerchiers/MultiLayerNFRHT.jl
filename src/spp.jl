# File for computing Surface-Plasmon-Polariton dispersion relation

" Parallel wavevector for SPP mode on a single interface "
kspp(eps1,eps2,w) = w/c0*sqrt(eps1*eps2/(eps1+eps2))


" Dispersion relation for SPP modes on a film "

function kk(kpar,eps,w)
    return sqrt(kpar^2 - eps*(w/c0)^2)
end

function kspp(kpar,a,eps1,eps2,eps3,w)
    k1 = kk(kpar,eps1,w)
    k2 = kk(kpar,eps2,w)
    k3 = kk(kpar,eps3,w)

    f1 = (k1/eps1 + k2/eps2)/(k1/eps1-k2/eps2)
    f2 = (k1/eps1 + k3/eps3)/(k1/eps1-k3/eps3)

    return exp(-4.0*k1*a)-f1*f2
end


function solve_kpar(a,eps1,eps2,eps3,w)
    f(kpar) = kssp(kpar,a,eps1,eps2,eps3,w)
    return fzero(f,w/c0,order = 16)
end

function dispersion_relation(a,material1 :: OptProp, material2 :: OptProp, material3 :: OptProp )
    wv = collect(linspace(1.0e13,1.0e15,1000))
    eps1 = permittivity.([material1],wv)
    eps2 = permittivity.([material2],wv)
    eps3 = permittivity.([material3],wv)
    kpar = zeros(Float64,1000)

    for i=1:1000
        kpar[i] = solve_kpar(a,eps1[i],eps2[i],eps3[i],wv[i])
    end

    return kpar
end
