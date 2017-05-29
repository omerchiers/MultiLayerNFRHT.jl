
using RadiativeHeat
using Plots, Cubature

function test1()
    wv = 1.0e15
    ang = collect(linspace(0.0,pi*0.5,1000))
    b = Bulk(Cst(),Cst(5.0+im*0.0))
    eps1 = permittivity(b.ep1,wv)
    kx = compute_kx.(ang,[eps1],[wv])
    RTE=zeros(Float64,1000)
    RTM=zeros(Float64,1000)

    for i=1:1000
        (rte,tte) = rt(b, te(), kx[i] ,wv)
        (rtm,ttm) = rt(b, tm(), kx[i] ,wv)
         RTE[i] = abs(rte)^2
         RTM[i] = abs(rtm)^2
    end
    plot(ang./pi.*180,RTE,xlim=(0,90),ylim=(0,1))
    plot!(ang./pi.*180,RTM,xlim=(0,90),ylim=(0,1))
end


function test2()
    wv = 1.0e15
    ang = collect(linspace(0.0,pi*0.5,1000))
    ml = [Tmp.Layer(Tmp.Cst()) ; Tmp.Layer(Tmp.Cst(5.0+im*3.0))]
    eps1 = Tmp.permittivity(ml[1].material,wv)
    kx = Tmp.compute_kx.(ang,[eps1],[wv])
    RTE=zeros(Float64,1000)
    RTM=zeros(Float64,1000)

    for i=1:1000
        (rte,tte) = Tmp.rt(ml, Tmp.te(), kx[i] ,wv)
        (rtm,ttm) = Tmp.rt(ml, Tmp.tm(), kx[i] ,wv)
         RTE[i] = abs(rte)^2
         RTM[i] = abs(rtm)^2
    end
    plot(ang./pi.*180,RTE,xlim=(0,90),ylim=(0,1))
    plot!(ang./pi.*180,RTM,xlim=(0,90),ylim=(0,1))
end

function test3()
    wv = 1.0e15
    ang = collect(linspace(0.0,pi*0.5,1000))
    ml = [Tmp.Layer(Tmp.Cst(5.0+im*3.0)) ; Tmp.Layer(Tmp.Cst(10.0+im*3.0))]
    eps1 = Tmp.permittivity(ml[1].material,wv)
    kx = Tmp.compute_kx.(ang,[eps1],[wv])
    RTE=zeros(Float64,1000)
    RTM=zeros(Float64,1000)

    for i=1:1000
        (rte,tte) = Tmp.rt(ml, Tmp.te(), kx[i] ,wv)
        (rtm,ttm) = Tmp.rt(ml, Tmp.tm(), kx[i] ,wv)
         RTE[i] = abs(rte)^2
         RTM[i] = abs(rtm)^2
    end
    plot(ang./pi.*180,RTE,xlim=(0,90),ylim=(0,1))
    plot!(ang./pi.*180,RTM,xlim=(0,90),ylim=(0,1))
end


function test4()
    wv = 1.0e13
    ang = collect(linspace(0.0,pi*0.5,1000))
    ml = [Tmp.Layer(Tmp.Cst()) ; Tmp.Layer(Tmp.Cst(5.0+im*3.0),1.0e-6) ; Tmp.Layer(Tmp.Cst(10.0+im*3.0))]
    eps1 = Tmp.permittivity(ml[1].material,wv)
    kx = Tmp.compute_kx.(ang,[eps1],[wv])
    RTE=zeros(Float64,1000)
    RTM=zeros(Float64,1000)

    for i=1:1000
        (rte,tte) = Tmp.rt(ml, Tmp.te(), kx[i] ,wv)
        (rtm,ttm) = Tmp.rt(ml, Tmp.tm(), kx[i] ,wv)
         RTE[i] = abs(rte)^2
         RTM[i] = abs(rtm)^2
      #   println("TTE=" ,abs(tte)^2,"TTM=" ,abs(ttm)^2)
    end
    plot(ang./pi.*180,RTE,xlim=(0,90),ylim=(0,1))
    plot!(ang./pi.*180,RTM,xlim=(0,90),ylim=(0,1))
end

function test5()
    wv = 1.0e13
    th = 1.0e-6
    ang = collect(linspace(0.0,pi*0.5,1000))
    ml = [Tmp.Layer(Tmp.Cst()) ; Tmp.Layer(Tmp.Cst(5.0+im*3.0),th) ; Tmp.Layer(Tmp.Cst(10.0+im*3.0))]
    eps1 = Tmp.permittivity(ml[1].material,wv)
    kx = Tmp.compute_kx.(ang,[eps1],[wv])
    RTE=zeros(Float64,1000)
    RTM=zeros(Float64,1000)

    for i=1:1000
        (rte01,tte) = Tmp.rt(Tmp.Bulk(Tmp.Cst(),Tmp.Cst(5.0+im*3.0)), Tmp.te(), kx[i] ,wv)
        (rtm01,ttm) = Tmp.rt(Tmp.Bulk(Tmp.Cst(),Tmp.Cst(5.0+im*3.0)), Tmp.tm(), kx[i] ,wv)

        (rte12,tte) = Tmp.rt(Tmp.Bulk(Tmp.Cst(5.0+im*3.0),Tmp.Cst(10.0+im*3.0)), Tmp.te(), kx[i] ,wv)
        (rtm12,ttm) = Tmp.rt(Tmp.Bulk(Tmp.Cst(5.0+im*3.0),Tmp.Cst(10.0+im*3.0)), Tmp.tm(), kx[i] ,wv)

        kz = Tmp.compute_kz(kx[i],5.0+im*3.0,wv)

        rte = (rte01 +rte12*exp(-im*2.0*kz*th))/(1.0+rte01*rte12*exp(-im*2.0*kz*th))
        rtm = (rtm01 +rtm12*exp(-im*2.0*kz*th))/(1.0+rtm01*rtm12*exp(-im*2.0*kz*th))

        RTE[i] = abs(rte)^2
        RTM[i] = abs(rtm)^2
      #   println("TTE=" ,abs(tte)^2,"TTM=" ,abs(ttm)^2)
    end
    plot(ang./pi.*180,RTE,xlim=(0,90),ylim=(0,1))
    plot!(ang./pi.*180,RTM,xlim=(0,90),ylim=(0,1))
end


function test6(ang)
    wv = collect(linspace(1.0e13,1.0e15,1000))
    ml = [Tmp.Layer(Tmp.Cst()) ; Tmp.Layer(Tmp.Cst(5.0+im*3.0),1.0e-6) ; Tmp.Layer(Tmp.Cst(10.0+im*3.0))]
    eps1 = Tmp.permittivity.([ml[1].material],wv)
    kx = Tmp.compute_kx.([ang],eps1,wv)
    RTE=zeros(Float64,1000)
    RTM=zeros(Float64,1000)

    for i=1:1000
        (rte,tte) = Tmp.rt(ml, Tmp.te(), kx[i] ,wv[i])
        (rtm,ttm) = Tmp.rt(ml, Tmp.tm(), kx[i] ,wv[i])
         RTE[i] = abs(rte)^2
         RTM[i] = abs(rtm)^2
    end
    plot(wv,RTE,xlim=(1e13,1e15),ylim=(0,1))
    plot!(wv,RTM,xlim=(1e13,1e15),ylim=(0,1))
    xaxis!(:log10)


end

function test7(T,th)

    wv = collect(linspace(1.0e13,1.0e15,1000))
    ml = [Tmp.Layer(Tmp.Cst()) ; Tmp.Layer(Tmp.Si(),th) ; Tmp.Layer(Tmp.Al())]
    em = zeros(Float64,1000)
    for i=1:1000
       em[i] = Tmp.emissivity_w(ml, wv[i])
    end
    plot(wv,em,xlim=(1e13,1e15),ylim=(0,0.14))
    xaxis!(:log10)

end

function test8(T1,T2,dist,dim)
    wv = collect(linspace(1.0e13,1.0e15,dim))
    ml = [Tmp.Layer(Tmp.Al())]
    gapl = Tmp.Layer(Tmp.Cst(),dist)
    qw         = zeros(Float64,dim)
    qw_ev_te   = zeros(Float64,dim)
    qw_ev_tm   = zeros(Float64,dim)
    qw_prop_te = zeros(Float64,dim)
    qw_prop_tm = zeros(Float64,dim)

    for i=1:dim
       qw_ev_te[i] = Tmp.heat_transfer_w(Tmp.Evanescent() ,ml, ml, gapl ,Tmp.te() ,wv[i], T1,T2)
       qw_ev_tm[i] = Tmp.heat_transfer_w(Tmp.Evanescent() ,ml, ml, gapl ,Tmp.tm() ,wv[i], T1,T2)
       qw_prop_te[i] = Tmp.heat_transfer_w(Tmp.Propagative() ,ml, ml, gapl ,Tmp.te() ,wv[i], T1,T2)
       qw_prop_tm[i] = Tmp.heat_transfer_w(Tmp.Propagative() ,ml, ml, gapl ,Tmp.tm() ,wv[i], T1,T2)
    end
    qw = qw_ev_te .+ qw_ev_tm .+ qw_prop_te .+qw_prop_tm

    plot(wv,qw,xlim=(1e13,1e15))
    plot!(wv,qw_ev_te,xlim=(1e13,1e15))
    plot!(wv,qw_ev_tm,xlim=(1e13,1e15))
    plot!(wv,qw_prop_te,xlim=(1e13,1e15))
    plot!(wv,qw_prop_tm,xlim=(1e13,1e15))

    xaxis!(:log10)
    yaxis!(:log10)

 end

 function test9(dist)
     wv = 1.0e15
     kx = collect(linspace(wv/Tmp.c0+1.0e-10,8.0*wv/Tmp.c0,1000))
     ml = [Tmp.Layer(Tmp.Al())]
     gapl  = Tmp.Layer(Tmp.Cst(),dist)
     tkx_w = zeros(Float64,1000)
     for i=1:1000
         tkx_w[i] = kx[i]*Tmp.transmission_kx_w(Tmp.Evanescent() ,ml,ml,gapl,Tmp.tm(),kx[i],wv)
         println(tkx_w[i])
    end
    plot(kx,tkx_w)
    xaxis!(:log10)
    yaxis!(:log10)
  end

  function test10(T1,T2,dim)
      ml = Tmp.Layer(Tmp.Al())
      dist = collect(logspace(-8,-4,dim))
      q = zeros(Float64,dim,5)

      dt = now()
      file_name = "data/"string(dt)*"_Al_Al_q_vs_dist.dat"
      f = open(file_name, "a")

      for i=1:dim
          #println("gap separation = ",dist[i])
          gapl   = Tmp.Layer(Tmp.Cst(),dist[i])
          q[i,:] = Tmp.total_heat_transfer(ml, ml, gapl, T1,T2)
          output_vec = [dist[i] q[i,:]']
          write(f, "$(output_vec) \r\n")
     end
    close(f)


     plot(dist,q)
     xaxis!(:log10)
     yaxis!(:log10)
   end

   function test11(dim)

       dist = collect(logspace(-8,-4,dim))
       q = zeros(Float64,dim,5)

       dt = now()
       file_name = "data/"string(dt)*"_Al_Al_q_vs_dist.dat"
       f = open(file_name, "a")
           for i=1:dim
               println("gap separation = ",dist[i])
               q[i,:] = rand(1,5)
               output_vec = [dist[i] q[i,:]']
               write(f, "$(output_vec) \r\n")
               close(f)
           end

      plot(dist,q)
      xaxis!(:log10)
      yaxis!(:log10)
    end

" compute kx-omega map of emissivity"
    function test12(th,dimw,dimk)
      ml = [Tmp.Layer(Tmp.Cst()) ; Tmp.Layer(Tmp.Sic(),th) ; Tmp.Layer(Tmp.Al())]
      wv = collect(linspace(1.0e13,1.0e15,dimw))
      em = zeros(Float64,dimw*dimk,3)
      count = 1
      for i = 1:dimw
        kx = collect(linspace(1.0e4,wv[i]/Tmp.c0,dimk))
        for j = 1:dimk
            em[count,1] = kx[j]
            em[count,2] = wv[i]
            em[count,3] = Tmp.emissivity_kx_w(ml, kx[j], wv[i])
            count += 1
        end
      end
      output_file = "data/emissivity/omega_kx_th=1000nm_SiC_Al.dat"
      writedlm(output_file,em)

    end

" compute dispersion relation SPP for Air-SiC interface"
    function test13(dim)
          wv = collect(linspace(1.0e13,1.0e15,dim))
          eps1 = Tmp.permittivity.([Tmp.Cst()],wv)
          eps2 = Tmp.permittivity.([Tmp.Sic()],wv)
          kpar = Tmp.kssp.(eps1,eps2,wv)
          kl   = real(sqrt.(eps2)).*wv/Tmp.c0
          output_data = [real(kpar) imag(kpar) wv kl ]
          output_file = "data/SPP/Air_SiC.dat"
          writedlm(output_file,output_data)
    end

    " compute dispersion relation SPP for SiC-Al interface"
    function test14(dim)
          wv = collect(linspace(1.0e13,1.0e15,dim))
          eps1 = Tmp.permittivity.([Tmp.Sic()],wv)
          eps2 = Tmp.permittivity.([Tmp.Al()],wv)
          kpar = Tmp.kssp.(eps1,eps2,wv)
          kl   = real(sqrt.(eps2)).*wv/Tmp.c0
          output_data = [real(kpar) imag(kpar) wv kl]
          output_file = "data/SPP/SiC_Al.dat"
          writedlm(output_file,output_data)
    end

" compute monocromatic hemispherical emissivity"
    function test15(dim,th)
          ml = [Tmp.Layer(Tmp.Cst()) ; Tmp.Layer(Tmp.Sic(),th) ; Tmp.Layer(Tmp.Al())]
          wv = collect(linspace(1.0e13,1.0e15,dim))
          em = Tmp.emissivity_w.([ml], wv)
          output_data = [wv em]
          output_file = "data/emissivity/SiC_th=1000nm_Al.dat"
          writedlm(output_file,output_data)
    end

" compute total emissivity and compare with classical formula for far-field transfer"
    function test16(th,T1,T2)
      ml = [Tmp.Layer(Tmp.Cst()) ; Tmp.Layer(Tmp.Sic(),th) ; Tmp.Layer(Tmp.Al())]
      em_ml = Tmp.emissivity(ml, T1)

      b = Tmp.Bulk(Tmp.Cst(),Tmp.Al())
      em_b = Tmp.emissivity(b, T2)
      hf = Tmp.farfield_transfer(em_ml , em_b, T1, T2)

      println("emissivity multilayer = ",em_ml)
      println("emissivity aluminium = ",em_b)
      println("heat flux = ",hf ,"W/m^2")

    end

 " compute total emissivity for two bulk media and compare with classical formula for far-field transfer"
        function test17(T1,T2)
          b1 = Tmp.Bulk(Tmp.Cst(),Tmp.Sic())
          em_ml = Tmp.emissivity(b1, T1)

          b2 = Tmp.Bulk(Tmp.Cst(),Tmp.Al())
          em_b = Tmp.emissivity(b2, T2)
          hf = Tmp.farfield_transfer(em_ml , em_b, T1, T2)

          println("emissivity SiC = ",em_ml)
          println("emissivity aluminium = ",em_b)
          println("heat flux = ",hf ,"W/m^2")

        end

" heat_flux with anonymous function method"
    function test18(tol)
        b1  = Layer(Al())
        gap = Layer(Cst(),1.0e-8)
        #@time heat_flux_integrand(Evanescent(),b1,b1,gap,te(),300.0,tol)
        @time q1 = heat_flux(Evanescent(),b1,b1,gap,te(),300.0; tol1=1e-8 , tol2=1e-8)
        println(q1)
        @time q2 = heat_flux2(Evanescent(),b1,b1,gap,te(),300.0, 1e-5)
        println(q2)
    end

" heat_transfer with anonymous function method"
    function test19(tol)
        b1  = Layer(Al())
        gap = Layer(Cst(),1.0e-8)
        #@time heat_flux_integrand(Evanescent(),b1,b1,gap,te(),300.0,tol)
        @time q1 = heat_transfer(Evanescent(),b1,b1,gap,te(),400.0 , 300.0  ; tol1=tol,tol2=tol)
        println("heat_transfer  =",q1)
        @time q2 = heat_transfer2(Evanescent(),b1,b1,gap,te(),400.0 , 300.0 ; tol1=tol,tol2=tol)
        println("heat_transfer2 =",q2)
    end

 " total heat_transfer as a function of separation distance"
        function test20(tl1,tl2)
            b1  = Layer(Al())
            dist = logspace(-8,-4,40)
            q = zeros(Float64,40,5)
            fid = open("q_vs_dist.dat", "a")
            #@time heat_flux_integrand(Evanescent(),b1,b1,gap,te(),300.0,tol)
            for i=40:40
              gap = Layer(Cst(),dist[i])
              @time q[i,1:5] = total_heat_transfer(b1,b1,gap, 400.0 , 300.0  , tl1,tl2)
              println("distance in m  =",dist[i])
              println("heat_transfer  =",q[i,1:5])
              writedlm(fid,[dist[i], q[i,1:5]'])
            end
            close(fid)
        end

  " heat_transfer_w"
  function test21(tol)
     wv  = collect(logspace(13,15,1000))
     b1  = Layer(Al(),0.0)
     gap = Layer(Cst(),1.0e-4)
     #@time heat_flux_integrand(Evanescent(),b1,b1,gap,te(),300.0,tol)
     q = zeros(Float64,1000)
     qp = zeros(Float64,1000)
     for i=1:1000
        for f in [ Evanescent(), Propagative()]
            for p in [te(),tm()]
                q[i] += heat_transfer_w(f ,b1 ,b1, gap ,p ,wv[i], 400.0,300.0;toler=tol)
             end
        end
        qp[i] = planck(wv[i],400.0)-planck(wv[i],300.0)
      end
      qtot=trapz(wv,q)
      println(qtot)
      plot(wv,q,xaxis=:log10,yaxis=:log10)
      plot!(wv,qp,xaxis=:log10,yaxis=:log10)

  end

  " transmission_w"
  function test22(tol)
     wv  = collect(linspace(1e13,1e15,1000))
     b1  = Layer(Cst(1.00001+im*0.001))
     gap = Layer(Cst(),1.0e-2)
     #@time heat_flux_integrand(Evanescent(),b1,b1,gap,te(),300.0,tol)
     q = zeros(Float64,1000)
     for i=1:1000
        for f in [Evanescent(), Propagative()]
            for p in [te(),tm()]
                q[i] += transmission_w(f ,b1 ,b1, gap ,p ,wv[i],tol)
             end
        end
      end
      println(q[1000])
      plot(wv,q,xaxis=:log10,yaxis=:log10)
    #  plot!(wv,wv.^2/pi/c0,xaxis=:log10,yaxis=:log10)

  end

  " total flux black body sigma*T^4 "
  function test23(T)
     wmin = wien(T)*0.2
     wmax = wien(T)*20.0
     println("wmin =",wmin)
     println("wmax =",wmax)
     wv  = collect(linspace(wmin,wmax,1000))

     qp  = zeros(Float64,1000)
     qbe = zeros(Float64,1000)

     for i=1:1000
           qbe[i] = bose_einstein(wv[i],T)*wv[i]^2/c0^2/4.0/pi^2
           qp[i]  = planck(wv[i],T)
     end
      qtotbe = trapz(wv,qbe)
      qtotp  = trapz(wv,qp)
      println(qtotbe/sigma/T^4)
      println(qtotp/sigma/T^4)
      plot(wv,qp,xaxis=:log10,yaxis=:log10)
      plot!(wv,qbe,xaxis=:log10,yaxis=:log10)

  end

  " Check if transmission_kx_w gives 2 for a black body"
  function test24()
     wv  = collect(linspace(1e13,1e15,1000))
     b1  = Layer(Cst(1.00001+im*0.001))
     gap = Layer(Cst(),1.0e-5)
     #@time heat_flux_integrand(Evanescent(),b1,b1,gap,te(),300.0,tol)
     q = zeros(Float64,1000)
     for i=1:1000
        for f in [Evanescent(), Propagative()]
            for p in [te(),tm()]
                q[i] += transmission_kx_w(f,b1,b1,gap,p,1e4,wv[i])
             end
        end
      end
      println(q[1000])
      plot(wv,q,xaxis=:log10)
    #  plot!(wv,wv.^2/pi/c0,xaxis=:log10,yaxis=:log10)

  end

  " Check convergence of evanescent part"

  function test25(tol,dist)
     wv  = collect(linspace(1e13,1e15,1000))
     b1  = Layer(Al())
     gap = Layer(Cst(),dist)
     #@time heat_flux_integrand(Evanescent(),b1,b1,gap,te(),300.0,tol)
     q = zeros(Float64,1000)
     q2 = zeros(Float64,1000)

     for i=1:1000
        for f in [Evanescent()]
            for p in [te(),tm()]
                q[i] += transmission_w(f ,b1 ,b1, gap ,p ,wv[i],tol)
                q2[i] += transmission_w2(f ,b1 ,b1, gap ,p ,wv[i],tol)
             end
        end
      end
      println("q = ",q[1000])
      println("q2 = ",q2[1000])
      plot(wv,q,xaxis=:log10,yaxis=:log10)
      plot!(wv,q2,xaxis=:log10,yaxis=:log10)

    #  plot!(wv,wv.^2/pi/c0,xaxis=:log10,yaxis=:log10)

  end


  " Check exponential factor in evanescent contribution"
  function test26(dist,w)
      #dist  = collect(logspace(-9,-6,100))
      kxv   = collect(linspace(w/c0,1e4*w/c0,1000))

      #@time heat_flux_integrand(Evanescent(),b1,b1,gap,te(),300.0,tol)
      exp_val = zeros(Float64,1000)
      fac_val = zeros(Float64,1000)
      b1  = Layer(Al())
      gap = Layer(Cst(),dist)
      for i=1:1000
          (r_21,t)= rt([Layer(Cst());b1],te(),kxv[i],w)
          r_23    = r_21
          k2z     = compute_kz(kxv[i],permittivity(Cst(),w),w)
          exp_val[i] = exp(-2.0*imag(k2z)*gap.thickness)
          exp_val2   = exp(2.0*im*k2z*gap.thickness)            :: Complex{Float64}
          fac_val[i] = imag(r_21)*imag(r_23)/abs(1.0-r_21*r_23*exp_val2)^2
      end
      #scatter(dist,exp_val,xaxis=:log10)

      plot(kxv,exp_val.*fac_val,yaxis=:log10)
      plot!(kxv,exp_val.*fac_val.*kxv,yaxis=:log10)
      plot!(kxv,1./kxv.^2,yaxis=:log10)
    #  plot!(wv,wv.^2/pi/c0,xaxis=:log10,yaxis=:log10)
  end

  " total heat_transfer as a function of separation distance with integration version"
  function test27(tol)
     wv  = collect(linspace(1e13,1e15,1000))
     dist= collect(logspace(-8,-4,100))
     b1  = Layer(Al(),0.0)
     b2  = [Layer(Sic(),600.0e-9) ; Layer(Al(),0.0)]
     #@time heat_flux_integrand(Evanescent(),b1,b1,gap,te(),300.0,tol)
     qtot = zeros(Float64,100)
     for i=1:100
        q = zeros(Float64,1000)
        gap = Layer(Cst(),dist[i])
        for j=1:1000
          for f in [ Evanescent(), Propagative()]
              for p in [te(),tm()]
                  q[j] += heat_transfer_w(f ,b1 ,b2, gap ,p ,wv[j], 400.0, 300.0 ;toler=tol)
               end
          end
        end
      qtot[i]=trapz(wv,q)
      println("distance =",dist[i])
      println("qtot =",qtot[i])
      end

      plot(dist,qtot,xaxis=:log10,yaxis=:log10)

  end

  " test for integration over frequency from 0 to infinity"
  function test28(tolw,tolkx)
     dist= collect(logspace(-8,-4,100))
     b1  = Layer(Al(),0.0)
     #@time heat_flux_integrand(Evanescent(),b1,b1,gap,te(),300.0,tol)
     qtot = zeros(Float64,100)
     for i=1:100
        valt = 0.0
        gap = Layer(Cst(),dist[i])
        for f in [ Evanescent(), Propagative()]
            for p in [te(),tm()]
                ht(w) = heat_transfer_w(f ,b1 ,b1, gap ,p ,w, 400.0,300.0;toler=tolkx)
                ht2(u) = ht(u*kb/ħ)
                ht3(t) = ht2(t/(1.0-t))/(1.0-t)^2
                (val,err) = hquadrature(ht3, 0.0 , 1.0 ; reltol=tolw, abstol=0, maxevals=0)
                valt +=val
             end
        end
      qtot[i]=valt*kb/ħ
      println("distance =",dist[i])
      println("qtot =",qtot[i])
      end

      plot(dist,qtot,xaxis=:log10,yaxis=:log10)

  end


  " total heat_transfer as a function of separation distance with trapezium integration version between Al-Diel"
  function test29(tol)
     file_name = "data/heat_transfer/Al-SiC/Al_SiC_q_vs_dist.dat"
     fn = open(file_name, "a")
     wv  = collect(linspace(1e12,1e15,2000))
     dist= collect(logspace(-8,-4,100))
     b1  = Layer(Al(),0.0)
     b2  = Layer(Al(),0.0)
     #@time heat_flux_integrand(Evanescent(),b1,b1,gap,te(),300.0,tol)
     qtot = zeros(Float64,100)
     for i=1:length(dist)
        q = zeros(Float64,length(wv))
        gap = Layer(Cst(),dist[i])
        for j=1:length(wv)
          for f in [ Evanescent(), Propagative()]
              for p in [te(),tm()]
                  q[j] += heat_transfer_w(f ,b1 ,b2, gap ,p ,wv[j], 400.0, 300.0 ;toler=tol)
               end
          end
        end
      qtot[i]=trapz(wv,q)
      output_vec = [dist[i] qtot[i]]
      write(fn, "$(output_vec) \r\n")
      println("distance =",dist[i])
      println("qtot =",qtot[i])
      end
      close(fn)
      plot(dist,qtot,xaxis=:log10,yaxis=:log10)

  end
