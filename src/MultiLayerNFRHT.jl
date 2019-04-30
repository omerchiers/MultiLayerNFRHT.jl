__precompile__(false)

module MultiLayerNFRHT

# Dependencies
using QuadGK, Roots, Reexport,HCubature
    @reexport using MyPhysicalConstants
    @reexport using OpticalProperties


# export functions
export bose_einstein,wien,lambda_wien,farfield_transfer,
       emissivity_kx_w,emissivity_kx,emissivity_w,emissivity,emissivity_fraction,
       transmission_kx_w,transmission_w,total_transmission_kx_w,
       total_transmission_map,
       compute_kz,compute_kx,
       heat_flux,heat_flux2,heat_flux_w,heat_flux_integrand,
       heat_flux_integrand,heat_transfer2,
       heat_transfer,heat_transfer_w,
       total_heat_transfer,total_heat_transfer_double,total_heat_transfer_w,
       rt,planck,planck_fraction,
       kspp,
       unitconv,
       trapz

# export types
export Structure,Bulk,Layer,MultiLayer,
       TotalField,Evanescent,Propagative,
       Polarization,te,tm,
       BulkOrMultiLayer,LayerOrMultiLayer


include("rt_coefficients.jl")
include("emissivity.jl")
include("heat_transfer.jl")
include("spp.jl")
include("utils.jl")

end # module
