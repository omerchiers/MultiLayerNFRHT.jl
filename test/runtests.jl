using MultiLayerNFRHT
using Base.Test

# Effective medium models
bm = Bruggeman(Cst(3.0 + 0.5*im),Cst(5.0 + 1.5*im),0.0)
mg = MaxwellGarnett(Cst(3.0 + 0.5*im),Cst(5.0 + 1.5*im),0.0)

@test permittivity(bm,1e15) == permittivity(Cst(5.0 + 1.5*im),1e15)
@test permittivity(mg,1e15) == permittivity(Cst(5.0 + 1.5*im),1e15)

# Far-field black body heat transfer
bb_bulk = Layer(Cst(1.0 + 1e-5*im))
gap     = Layer(Cst(),1e-1)

@test total_heat_transfer(b,b,gap,1.0,0.0,1e8,1e13) ≈ [6.17311e-14, 6.17343e-14 , 2.83518e-8 , 2.83518e-8 , 5.67038e-8 ]
