# This file contains the functions for heat transfer between plane geometries


" Far-field heat transfer between two semi-infinite media,classical approximation "
function farfield_transfer(em1 , em2, T1, T2)

    return sigma*(T2^4-T1^4)/(1.0/em1 + 1.0/em2 - 1.0)

end

" Evanescent contribution to radiative heat transfer "
function transmission_kx_w(:: Evanescent ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization , kx ,w)

    b1 = prepend!([b1],[gap])
    b2 = prepend!([b2],[gap])

    (r_21,t)  = rt(b1,pol,kx,w) :: Tuple{Complex{Float64},Complex{Float64}}
    (r_23,t)  = rt(b2,pol,kx,w) :: Tuple{Complex{Float64},Complex{Float64}}

    k2z = compute_kz(kx,permittivity(gap.material,w),w) :: Complex{Float64}
    exp_val1 = exp(2.0*imag(k2z)*gap.thickness)         :: Float64
    exp_val2 = exp(2.0*im*k2z*gap.thickness)            :: Complex{Float64}

    return 4.0*exp_val1*imag(r_21)*imag(r_23)/abs(1.0-r_21*r_23*exp_val2)^2

end


" Propagative contribution to radiative heat transfer "
function transmission_kx_w(:: Propagative, b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization , kx ,w)

    b1 = prepend!([b1],[gap])
    b2 = prepend!([b2],[gap])

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
    t3_kx_w(v)  = t_kx_w(1.0 + v/(1.0-v))/(1.0-v)^2

    val :: Float64  = 0.0
    err :: Float64  = 0.0
   (val,err) = hquadrature(t3_kx_w, 0.0 , 1.0 ; reltol=tol, abstol=0, maxevals=0)
    return val
end


function transmission_w(field :: Propagative ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization ,w,tol)
    t_kx_w(kx)  = kx*transmission_kx_w(field ,b1,b2,gap,pol,kx,w)

    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = hquadrature(t_kx_w, 0.0, w/c0 ; reltol=tol, abstol=0, maxevals=0)
    return val
end

" Monocromatic heat flux "
function heat_flux_w(field :: TotalField ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization ,w, T)
    t_w  = transmission_w(field,b1,b2,gap,pol,w)
    be_w = bose_einstein(w,T)
    return be_w*t_w/4.0/pi^2
end

" Heat flux "
function heat_flux(field :: TotalField ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization , T; tol1=1e-8 , tol2=1e-8)
    q_w(w)  = transmission_w(field,b1,b2,gap,pol,w,tol1)
    q2_w(u) = u*q_w(u*kb*T/ħ)/(exp(u)-1.0)
    q3_w(t) = q2_w(t/(1.0-t))/(1.0-t)^2

    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = hquadrature(q3_w, 0.0 , 1.0 ; reltol=tol2, abstol=0, maxevals=0)
    return val*(kb*T)^2/ħ/4.0/pi^2
end

function heat_flux_integrand(field :: TotalField ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization , T; tol=1e-8)
    q_w(w)  = transmission_w(field,b1,b2,gap,pol,w,tol)
    q2_w(u) = u*q_w(u*kb*T/ħ)/(exp(u)-1.0)
    q3_w(t) = q2_w(t/(1.0-t))/(1.0-t)^2

    return q3_w(t)
end



function heat_flux2(field :: Propagative ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization , T,tol)

    function t_kx_w(u,v)
       # frequency
       w  = u*kb*T/ħ
       # kx
       kx = v
       if kx <= w/c0
           return kx*transmission_kx_w(field ,b1,b2,gap,pol,kx,w)
       elseif kx > w/c0
           return 0.0
       end
    end

    be_w(u) = u/(exp(u)-1.0)

    integr1(u,v) = t_kx_w(u,v)*be_w(u)
    integr2(x)   = integr1(x[1]/(1.0-x[1]),x[2]/(1.0-x[2]))/(1.0-x[1])^2/(1.0-x[2])^2

    #limits
      xmin = [0.0 ; 0.0]
      xmax = [1.0 ; 1.0]
      val :: Float64  = 0.0
      err :: Float64  = 0.0
      (val,err) = hcubature(integr2, xmin , xmax ; reltol=tol, abstol=0, maxevals=0)
    return  val*(kb*T)^2/ħ/4.0/pi^2
end

function heat_flux2(field :: Evanescent ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization , T,tol)
  function t_kx_w(u,v)
     # frequency
     w  = u*kb*T/ħ
     # kx
     kx = v
     if kx >= w/c0
         return kx*transmission_kx_w(field ,b1,b2,gap,pol,kx,w)
     elseif kx < w/c0
         return 0.0
     end
  end

  be_w(u) = u/(exp(u)-1.0)

  integr1(u,v) = t_kx_w(u,v)*be_w(u)
  integr2(x)   = integr1(x[1]/(1.0-x[1]),x[1]/(1.0-x[1])+x[2]/(1.0-x[2]))/(1.0-x[1])^2/(1.0-x[2])^2

  #limits
    xmin = [0.0 ; 0.0]
    xmax = [1.0 ; 1.0]
    val :: Float64  = 0.0
    err :: Float64  = 0.0
    (val,err) = hcubature(integr2, xmin , xmax ; reltol=tol, abstol=0, maxevals=0)
  return  val*(kb*T)^2/ħ/4.0/pi^2
end




" Monocromatic net heat transfer "
function heat_transfer_w(field :: TotalField ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization ,w, T1,T2)
    return heat_flux_w(field,b1,b2,gap,pol,w,T2) - heat_flux_w(field,b1,b2,gap,pol,w,T1)
end

" Net heat transfer "
function heat_transfer(field :: TotalField ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization, T1,T2)
    return heat_flux(field,b1,b2,gap,pol,T2) - heat_flux(field,b1,b2,gap,pol,T1)
end

function heat_transfer2(field :: TotalField ,b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer ,pol :: Polarization, T1,T2)
    return heat_flux(field,b1,b2,gap,pol,T2) - heat_flux(field,b1,b2,gap,pol,T1)
end


"Total net heat transfer"
function total_heat_transfer(b1 :: LayerOrMultiLayer, b2 :: LayerOrMultiLayer, gap :: Layer, T1,T2)
    q_evan_ev_te = heat_transfer(Evanescent(),b1,b2,gap,te(),T1,T2)
    q_evan_ev_tm = heat_transfer(Evanescent(),b1,b2,gap,tm(),T1,T2)

    q_prop_ev_te = heat_transfer(Propagative(),b1,b2,gap,te(),T1,T2)
    q_prop_ev_tm = heat_transfer(Propagative(),b1,b2,gap,tm(),T1,T2)

    q_tot = q_evan_ev_te+q_evan_ev_tm + q_prop_ev_te+q_prop_ev_tm
    q = [q_tot  q_evan_ev_te  q_evan_ev_tm  q_prop_ev_te q_prop_ev_tm]

    return q

end
