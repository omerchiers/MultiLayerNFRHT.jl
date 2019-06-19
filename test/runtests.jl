using MultiLayerNFRHT
using Test

# Structures
bb_bulk = Layer(Cst(1.0 + 1e-5*im))
gap     = Layer(Cst(),1e-1)
ml1     = [Layer(Cst()),Layer(Cst(1.0+1e-5*im))]
ml2     = [Layer(Cst()),Layer(Au())]
ml3     = [Layer(Cst()),Layer(Al())]
ml4     = [Layer(Cst()),Layer(SiC)]
ml5     = [Layer(Cst()),Layer(SiC,1e-7),Layer(Al())]
ml6     = [Layer(Cst()),Layer(SiC,1e-6),Layer(Al())]


@testset "kz wavevector component" begin
    perm = 5.0 + 0.0*im
    w    = 1e15
    kx   = w/c0
    @test imag(compute_kz(kx,1.0,w)) == 0.0
    @test imag(compute_kz(kx,perm,w)) == 0.0
end


@testset "reflection coefficient" begin

end


@testset "spectral directional emissivity" begin

end

@testset "spectral hemispherical emissivity" begin

end

@testset "total emissivity" begin
    @test emissivity(ml1,300.0) ≈ 0.9999999658493035   atol = 1e-16
    @test emissivity(ml2,300.0) ≈ 0.018976160729311516 atol = 1e-16
    @test emissivity(ml3,300.0) ≈ 0.012711414795024877 atol = 1e-16
    @test emissivity(ml4,300.0) ≈ 0.6103470572275264   atol = 1e-16
    @test emissivity(ml5,300.0) ≈ 0.016788049404917925  atol = 1e-16
    @test emissivity(ml6,300.0) ≈ 0.06548138619305238  atol = 1e-16
end

@testset "total heat transfer" begin
    @test total_heat_transfer(bb_bulk,bb_bulk,gap,1.0,0.0,1e8,1e13) ≈ [6.17311e-14, 6.17343e-14 , 2.83518e-8 , 2.83518e-8 , 5.67038e-8 ] atol=1e-13
end
