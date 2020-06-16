# This file contains the functions for heat transfer between plane geometries


" Far-field heat transfer between two semi-infinite media,classical approximation "
function farfield_transfer(em1 , em2, T1, T2)

    return sigma*(T2^4-T1^4)/(1.0/em1 + 1.0/em2 - 1.0)

end

"""
    incoherent_total_heat_transfer_w(b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer, T1,T2,w;tolkx=1e-6)

compute the total monocromatic radiative heat flux between planar surfaces in the far field,
considereing total incoherent radiation. This should only be considered correct (less than 5 % error) when the separation
distance between both objects is larger than 100*λ_wien.

#Arguments:
- `b1 :: LayerOrMultiLayer`: Multilayer description of body 1
- `b2 :: LayerOrMultiLayer`: Multilayer description of body 2
- `gap :: Layer`: Material and thickness of gap separation
- `T1`: temperature of body 1
- `T2`: temperature of body 2
- `w` : frequency in radHz
- `tolkx = 1e-6`: relative error on integration over kx
"""
function incoherent_total_heat_transfer_w(b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer, T1,T2,w;tolkx=1e-6)
   T1==0.0 ? l1=0.0 : l1 = lambda_wien(T1)
   T2==0.0 ? l2=0.0 : l2 = lambda_wien(T2)

   @assert gap.thickness >= 100*max(l1,l2)  "the gap size should at least be 100 times Wien's wavelength"

    em1(kx) = emissivity_kx_w([Layer(gap.material) ; b1], kx, w)
    em2(kx) = emissivity_kx_w([Layer(gap.material) ; b2], kx, w)

    tn(kx)  = kx*em1(kx)*em2(kx)*0.25
    td(kx)  = 1 - ( 1- em1(kx)*0.5)*( 1- em2(kx)*0.5)
    t(kx)   = tn(kx)/td(kx)
    t2(u)   = t(u*w/c0)

    (val,err) = quadgk(t2, 0,1; rtol=tolkx)
    return 2/(w/c0)*val*(planck(w,T1) - planck(w,T2))
end

"""
    incoherent_total_heat_transfer(b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer, T1, T2, w1, w2; tolkx=1e-6)

compute the total radiative heat flux between planar surfaces in the far field,
considereing total incoherent radiation. This should only be considered correct (less than 5 % error) when the separation
distance between both objects is larger than 100*λ_wien.

#Arguments:
- `b1 :: LayerOrMultiLayer`: Multilayer description of body 1
- `b2 :: LayerOrMultiLayer`: Multilayer description of body 2
- `gap :: Layer`: Material and thickness of gap separation
- `T1`: temperature of body 1
- `T2`: temperature of body 2
- `w` : frequency in radHz
- `tolkx = 1e-6`: relative error on integration over kx
- `tolw = 1e-6`: relative error on integration over the frequency w
"""
function incoherent_total_heat_transfer(b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer, T1,T2,w1,w2;tolkx=1e-6, tolw = 1e-6)
    u1 = w1*ħ/kb
    u2 = w2*ħ/kb
    ht(w) = incoherent_total_heat_transfer_w(b1 ,b2, gap , T1,T2, w ;tolkx=tolkx)
    ht2(u) = ht(u*kb/ħ)
    (val,err) = quadgk(ht2, u1 , u2 ; rtol=tolw)
    return val*kb/ħ
end


"""
    incoherent_total_heat_transfer(b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer, T1, T2; tolkx=1e-6)

compute integration over ω from 0 to +∞
"""
function incoherent_total_heat_transfer(b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer, T1,T2;tolkx=1e-6, tolw = 1e-6)
    ht(w) = incoherent_total_heat_transfer_w(b1 ,b2, gap , T1,T2, w ;tolkx=tolkx)
    ht2(u) = ht(u*kb/ħ)
    ht3(t) = ht2(t/(1-t))/(1-t)^2
    (val,err) = quadgk(ht3, 0 , 1 ; rtol=tolw)
    return val*kb/ħ
end



" Evanescent contribution to radiative heat transfer "
function transmission_kx_w(:: Evanescent ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization , kx ,w)

    gap2 = Layer(gap.material,0.0)
    b1   = prepend!([b1;],[gap2])
    b2   = prepend!([b2;],[gap2])

    (r_21,t)  = rt(b1,pol,kx,w) :: Tuple{Complex{Float64},Complex{Float64}}
    (r_23,t)  = rt(b2,pol,kx,w) :: Tuple{Complex{Float64},Complex{Float64}}

    k2z = compute_kz(kx,permittivity(gap.material,w),w) :: Complex{Float64}
    exp_val1 = exp(-2.0*imag(k2z)*gap.thickness)         :: Float64
    exp_val2 = exp(2.0*im*k2z*gap.thickness)            :: Complex{Float64}

    return 4.0*exp_val1*imag(r_21)*imag(r_23)/abs(1.0-r_21*r_23*exp_val2)^2

end


" Propagative contribution to radiative heat transfer "
function transmission_kx_w(:: Propagative, b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization , kx ,w)

    gap2 = Layer(gap.material,0.0)
    b1 = prepend!([b1;],[gap2])
    b2 = prepend!([b2;],[gap2])

    T_21      = 0.0
    T_23      = 0.0
    (r_21,t_21)  = rt(b1,pol,kx,w) :: Tuple{Complex{Float64},Complex{Float64}}
    (r_23,t_23)  = rt(b2,pol,kx,w) :: Tuple{Complex{Float64},Complex{Float64}}

    k2z = compute_kz(kx,permittivity(gap.material,w),w)  :: Complex{Float64}
    exp_val2 = exp(2.0*im*k2z*gap.thickness)             :: Complex{Float64}

    if imag(permittivity(substrate(b1),w)) == 0.0
        eps1,k1z,eps2,k2z = outer_media(b1,kx,w)
        T_21 = power_t(pol,eps1,eps2, k1z,k2z,t_21)
    end
    if imag(permittivity(substrate(b2),w)) == 0.0
        eps1,k1z,eps2,k2z = outer_media(b2,kx,w)
        T_23 = power_t(pol,eps1,eps2, k1z,k2z,t_23)
    end

    return (1.0 - abs(r_21)^2 - T_21)*(1.0 - abs(r_23)^2 - T_23)/abs(1.0-r_21*r_23*exp_val2)^2

end



" Monocromatic transmission "
function transmission_w(field :: Evanescent ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization ,w,tol)

    t_kx_w(kx)  = kx*transmission_kx_w(field ,b1,b2,gap,pol,kx,w)
    t2_kx_w(u)  = t_kx_w(u*w/c0)
    t3_kx_w(v)  = t2_kx_w(1 + v/(1-v))/(1-v)^2

   (val,err) = quadgk(t3_kx_w, 0,1; rtol=tol)
    return val*(w/c0)
end

function transmission_w(field :: Propagative ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization ,w,tol)
    t_kx_w(kx)  = kx*transmission_kx_w(field ,b1,b2,gap,pol,kx,w)
    t2_kx_w(u)  = t_kx_w(u*w/c0)

    (val,err) = quadgk(t2_kx_w, 0,1; rtol=tol)
    return val*(w/c0)
end

" Monocromatic heat flux "
function heat_flux_w(field :: TotalField ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization ,w, T;tol=1e-2)
    t_w  = transmission_w(field,b1,b2,gap,pol,w,tol)
    be_w = bose_einstein(w,T)
    return be_w*t_w/4.0/pi^2
end

" Heat flux "
function heat_flux(field :: TotalField ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization , T; w1=1e13,w2=1e15, tol1=1e-8 , tol2=1e-8)
    q_w(w)  = heat_flux_w(field ,b1 , b2 , gap  ,pol , w ,T ; tol=tol1)
    q1_w(u) = q_w(u*kb*T/ħ)
    u1=w1*ħ/kb/T
    u2=w2*ħ/kb/T

    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = quadgk(q1_w, u1 , u2 ; rtol=tol2)

    return val*(kb*T)/ħ
end




" Monocromatic net heat transfer "
function heat_transfer_w(field :: TotalField ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization ,w, T1,T2;toler=1e-2)
    return (bose_einstein(w,T1)-bose_einstein(w,T2))*transmission_w(field,b1,b2,gap,pol,w,toler)/4.0/pi^2
end

" Net heat transfer "
function heat_transfer(field :: TotalField ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization, T1,T2 ; tol1=1e-8,tol2=1e-8)
    return heat_flux(field,b1,b2,gap,pol,T2 ; tol1=tol1 , tol2=tol2) - heat_flux(field,b1,b2,gap,pol,T1 ; tol1=tol1 , tol2=tol2)
end



function total_heat_transfer_w(b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer, T1,T2,w;tolkx=1e-6,tolw=1e-6)
    valt = 0.0
    cnt  = 0
    q    = zeros(Float64,5)
    for f in ( Evanescent(), Propagative())
        for p in (te(),tm())
            cnt += 1
            valt  += heat_transfer_w(f ,b1 ,b2, gap ,p ,w, T1,T2;toler=tolkx)
         end
    end
    return valt
end

"Total net heat transfer"
function total_heat_transfer(b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer, T1,T2;tolkx=1e-6,tolw=1e-6)
    valt = 0.0
    cnt  = 0
    q    = zeros(Float64,5)
    for f in ( Evanescent(), Propagative())
        for p in (te(),tm())
            cnt  += 1
            ht(w) = heat_transfer_w(f ,b1 ,b2, gap ,p ,w, T1,T2;toler=tolkx)
            ht2(u) = ht(u*kb/ħ)
            ht3(t) = ht2(t/(1-t))/(1-t)^2
            (val,err) = quadgk(ht3, 0 , 1 , rtol=tolw)
            valt  += val
            q[cnt] = val*kb/ħ
         end
    end
    q[5]=valt*kb/ħ

    return q
end

function total_heat_transfer(b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer, T1,T2,w1,w2;tolkx=1e-6,tolw=1e-6)
    valt = 0.0
    cnt  = 0
    q    = zeros(Float64,5)
    u1 = w1*ħ/kb
    u2 = w2*ħ/kb
    for f in (Evanescent(),Propagative())
        for p in (te(),tm())
            cnt  += 1
            ht(w) = heat_transfer_w(f ,b1 ,b2, gap ,p ,w, T1,T2;toler=tolkx)
            ht2(u) = ht(u*kb/ħ)
            (val,err) = quadgk(ht2, u1 , u2 ; rtol=tolw)
            valt  += val
            q[cnt] = val*kb/ħ
         end
    end
    q[5]=valt*kb/ħ

    return q
end

function total_heat_transfer_double(b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer, T1,T2,w1,w2;tol=1e-6)
    valt = 0.0
    cnt  = 0
    q    = zeros(Float64,5)
    u1 = w1*ħ/kb
    u2 = w2*ħ/kb
    for f in (Evanescent(), Propagative())
        for p in (te(),tm())
            cnt  += 1
            val = integrand_double(f, p, b1, b2, gap, T1,T2,w1,w2;tol=tol)
            valt  += val
            q[cnt] = val*(kb/ħ)^3/c0^2/4.0/pi^2
         end
    end
    q[5]=valt*(kb/ħ)^3/c0^2/4.0/pi^2

    return q
end

function integrand_double(field :: Propagative, pol :: Polarization, b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer, T1,T2,w1,w2;tol=1e-6)
    u1 = w1*ħ/kb
    u2 = w2*ħ/kb
    integr(x) = (bose_einstein(x[1]*kb/ħ,T1) - bose_einstein(x[1]*kb/ħ,T2))*x[2]*x[1]^2*transmission_kx_w(field ,b1,b2,gap,pol,x[2]*x[1]*kb/ħ/c0,x[1]*kb/ħ)
    (val,err) = hcubature(integr,[u1,0.0],[u2,1.0]; rtol=tol)
    return val
end

function integrand_double(field :: Evanescent, pol :: Polarization , b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer, T1,T2,w1,w2;tol=1e-6)
    u1 = w1*ħ/kb
    u2 = w2*ħ/kb
    integr(x)  = (bose_einstein(x[1]*kb/ħ,T1) - bose_einstein(x[1]*kb/ħ,T2))*x[2]*x[1]^2*transmission_kx_w(field ,b1,b2,gap,pol,x[2]*x[1]*kb/ħ/c0,x[1]*kb/ħ)
    integr2(y) = integr([y[1],1 + y[2]/(1-y[2])])/(1-y[2])^2
    (val,err)  = hcubature(integr2,[u1,0.0],[u2,1.0]; rtol=tol)
    return val
end


function total_transmission_kx_w(b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer, kx, w)
    valt = 0.0
    for p in (te(),tm())
        if kx <= w/c0
            valt += transmission_kx_w(Propagative(), b1, b2 , gap  ,p , kx ,w)
        else
            valt += transmission_kx_w(Evanescent(), b1, b2 , gap  ,p , kx ,w)
        end
    end
    return valt
end


function total_transmission_map(b1 :: LayerOrMultiLayer,
                                b2 :: LayerOrMultiLayer,
                                gap :: Layer,
                                kx :: AbstractArray, w :: AbstractArray)

    t = zeros(length(w),length(kx))
    for j = 1:length(kx)
        for i = 1:length(w)
            t[i,j] = total_transmission_kx_w(b1,b2,gap,kx[j],w[i])
        end
    end
    return t
end
