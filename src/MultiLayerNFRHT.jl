__precompile__(false)

module MultiLayerNFRHT

# Dependencies
using QuadGK, Roots, Reexport, HCubature
@reexport using MyPhysicalConstants
@reexport using OpticalProperties

# export functions
export bose_einstein,
    wien,
    lambda_wien,
    farfield_transfer,
    emissivity_kx_w,
    emissivity_kx,
    emissivity_w,
    emissivity,
    emissivity_fraction,
    transmission_kx_w,
    transmission_w,
    total_transmission_kx_w,
    total_transmission_map,
    compute_kz,
    compute_kx,
    heat_flux,
    heat_flux_w,
    heat_transfer,
    heat_transfer_w,
    total_heat_transfer,
    total_heat_transfer_double,
    total_heat_transfer_w,
    incoherent_total_heat_transfer,
    incoherent_total_heat_transfer_w,
    rt,
    planck,
    planck_fraction,
    kspp,
    unitconv,
    trapz

# export types
export Structure,
    Bulk,
    Layer,
    MultiLayer,
    Film,
    LayerOrMultiLayer,
    TotalField,
    Evanescent,
    Propagative,
    Polarization,
    te,
    tm

include("rt_coefficients.jl")
include("emissivity.jl")
include("heat_transfer.jl")
include("spp.jl")
include("utils.jl")

end # module
