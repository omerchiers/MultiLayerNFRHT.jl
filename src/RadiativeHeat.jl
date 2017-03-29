
module RadiativeHeat
using Plots

# export functions
export bose_einstein,wien,farfield_transfer,
       transmission_kx_w,transmission_w,
       heat_flux,heat_flux2,heat_flux_w,
       heat_flux_integrand,
       heat_transfer,heat_transfer_w,
       total_heat_transfer

#export types
export OptProp,Model,
       TotalField,Evanescent,Propagative,
       Polarization,te,tm,
       Structure,Bulk,Layer

# export constants
export ħ, kb, c0


include("physical_constants.jl")
include("optical_properties.jl")
include("rt_coefficients.jl")
include("emissivity.jl")
include("heat_transfer.jl")
include("SPP.jl")


end
