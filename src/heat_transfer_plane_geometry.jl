# This file contains the functions for heat transfer between plane geometries


" Far-field heat transfer between two semi-infinite media,classical approximation "
function farfield_transfer(em1 , em2, T1, T2)

    return sigma*(T2^4-T1^4)/(1.0/em1 + 1.0/em2 - 1.0)

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

    (r_21,t) = rt(b1,pol,kx,w) :: Tuple{Complex{Float64},Complex{Float64}}
    (r_23,t) = rt(b2,pol,kx,w) :: Tuple{Complex{Float64},Complex{Float64}}

    k2z = compute_kz(kx,permittivity(gap.material,w),w)  :: Complex{Float64}
    exp_val2 = exp(2.0*im*k2z*gap.thickness)             :: Complex{Float64}


    return (1.0-abs(r_21)^2)*(1.0-abs(r_23)^2)/abs(1.0-r_21*r_23*exp_val2)^2

end



" Monocromatic transmission "
function transmission_w(field :: Evanescent ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization ,w,tol)

    t_kx_w(kx)  = kx*transmission_kx_w(field ,b1,b2,gap,pol,kx,w)
    t2_kx_w(u)  = t_kx_w(u*w/c0)
    t3_kx_w(v)  = t2_kx_w(1.0 + v/(1.0-v))/(1.0-v)^2

    val :: Float64  = 0.0
    err :: Float64  = 0.0
   (val,err) = hquadrature(t3_kx_w, 0.0 , 1.0 ; reltol=tol, abstol=0, maxevals=0)
    return val*(w/c0)
end

function transmission_w(field :: Propagative ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization ,w,tol)
    t_kx_w(kx)  = kx*transmission_kx_w(field ,b1,b2,gap,pol,kx,w)
    t2_kx_w(u)  = t_kx_w(u*w/c0)

    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = hquadrature(t2_kx_w, 0.0 , 1.0 ; reltol=tol, abstol=0, maxevals=0)
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
    (val,err) = hquadrature(q1_w, u1 , u2 ; reltol=tol2, abstol=0, maxevals=0)

    return val*(kb*T)/ħ
end

function heat_flux_integrand(field :: TotalField ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization , T,tol)
    q_w(w)  = transmission_w(field,b1,b2,gap,pol,w,tol)
    q2_w(u) = u*q_w(u*kb*T/ħ)/(exp(u)-1.0)
    return q2_w
end



" Monocromatic net heat transfer "
function heat_transfer_w(field :: TotalField ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization ,w, T1,T2;toler=1e-2)
    return (bose_einstein(w,T1)-bose_einstein(w,T2))*transmission_w(field,b1,b2,gap,pol,w,toler)/4.0/pi^2
end

" Net heat transfer "
function heat_transfer(field :: TotalField ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization, T1,T2 ; tol1=1e-8,tol2=1e-8)
    return heat_flux(field,b1,b2,gap,pol,T2 ; tol1=tol1 , tol2=tol2) - heat_flux(field,b1,b2,gap,pol,T1 ; tol1=tol1 , tol2=tol2)
end

function heat_transfer2(field :: TotalField ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization, T1,T2 ; tol1=1e-8,tol2=1e-8)

    q1_w = heat_flux_integrand(field ,b1 , b2 , gap  ,pol , T1 , tol1)
    q2_w = heat_flux_integrand(field ,b1 , b2 , gap  ,pol , T2 , tol1)

    function q_w(t,v)
        v[1] = q1_w(t)
        v[2] = q2_w(t)
        return v
    end

    (val,err) = hquadrature(2,q_w, 0.0 , 1.0 ; reltol=tol2, abstol=1e-8, maxevals=0)
    return val[1]*(kb*T1)^2/ħ/4.0/pi^2-val[2]*(kb*T2)^2/ħ/4.0/pi^2

end


"Total net heat transfer"
function total_heat_transfer(b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer, T1,T2;tolkx=1e-6,tolw=1e-6)
    valt = 0.0
    cnt  = 0
    q    = zeros(Float64,5)
    for f in [ Evanescent(), Propagative()]
        for p in [te(),tm()]
            cnt  += 1
            ht(w) = heat_transfer_w(f ,b1 ,b2, gap ,p ,w, T1,T2;toler=tolkx)
            ht2(u) = ht(u*kb/ħ)
            ht3(t) = ht2(t/(1.0-t))/(1.0-t)^2
            (val,err) = hquadrature(ht3, 0.0 , 1.0 ; reltol=tolw, abstol=0, maxevals=0)
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
    for f in [ Evanescent(), Propagative()]
        for p in [te(),tm()]
            cnt  += 1
            ht(w) = heat_transfer_w(f ,b1 ,b2, gap ,p ,w, T1,T2;toler=tolkx)
            ht2(u) = ht(u*kb/ħ)
            (val,err) = hquadrature(ht2, u1 , u2 ; reltol=tolw, abstol=0, maxevals=0)
            valt  += val
            q[cnt] = val*kb/ħ
         end
    end
    q[5]=valt*kb/ħ

    return q
end

# total_heat_transfer(spectrum :: FrequencyRange,
#                     b1 :: LayerOrMultiLayer,
#                     b2 :: LayerOrMultiLayer,
#                     gap :: Layer,
#                     T1,T2;tolkx=1e-6,tolw=1e-6)
# = total_heat_transfer(b1 ,b2, gap ,T1,T2,0.1*wien(min(T1,T2)),10.0*wien(max(T1,T2));tolkx=1e-6,tolw=1e-6)
