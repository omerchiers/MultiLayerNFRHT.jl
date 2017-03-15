
using Tmp
using Plots

function test1()
    wv = 1.0e15
    ang = collect(linspace(0.0,pi*0.5,1000))
    b = Tmp.Bulk(Tmp.Cst(),Tmp.Cst(5.0+im*0.0))
    eps1 = Tmp.permittivity(b.ep1,wv)
    kx = Tmp.compute_kx.(ang,[eps1],[wv])
    RTE=zeros(Float64,1000)
    RTM=zeros(Float64,1000)

    for i=1:1000
        (rte,tte) = Tmp.rt(b, Tmp.te(), kx[i] ,wv)
        (rtm,ttm) = Tmp.rt(b, Tmp.tm(), kx[i] ,wv)
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
