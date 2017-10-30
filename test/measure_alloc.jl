using Tmp

function wrapper()
   b1  = Tmp.Layer(Tmp.Al(),0.0)
   gap = Tmp.Layer(Tmp.Cst(),1.0e-7)
   println(Tmp.heat_transfer(Tmp.Evanescent() ,b1, b1 , gap ,Tmp.te(), 300.0,400.0))
   Profile.clear_malloc_data()
   println(Tmp.heat_transfer(Tmp.Evanescent() ,b1, b1 , gap ,Tmp.te(), 300.0,400.0))
end

wrapper()
