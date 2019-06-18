using MultiLayerNFRHT
using Base.Test

# Far-field black body heat transfer
bb_bulk = Layer(Cst(1.0 + 1e-5*im))
gap     = Layer(Cst(),1e-1)

@test total_heat_transfer(b,b,gap,1.0,0.0,1e8,1e13) ≈ [6.17311e-14, 6.17343e-14 , 2.83518e-8 , 2.83518e-8 , 5.67038e-8 ]

# emissivity single bodies
ml1 = [Layer(Cst()),Layer(Cst(1.0+1e-5*im))]
@test emissivity(ml1,300.0) ≈ 0.9999999658493035

ml2 = [Layer(Cst()),Layer(Au())]
@test emissivity(ml2,300.0) ≈ 0.018976160729311516

ml3 = [Layer(Cst()),Layer(Al())]
@test emissivity(ml3,300.0) ≈ 0.012711414795024877

ml4 = [Layer(Cst()),Layer(Sic())]
@test emissivity(ml4,300.0) ≈ 0.6103470572275264
