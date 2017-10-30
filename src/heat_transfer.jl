
abstract type TotalField end
struct Evanescent  <: TotalField end
struct Propagative <: TotalField end

# abstract Spectrum
# immutable Frequency <: Spectrum end
# immutable Wavelength <: Spectrum end
#
# immutable FrequencyRange <: Spectrum
#     w1 :: Float64
#     w2 :: Float64
# end
#
# FrequencyRange() = FrequencyRange(0.0,0.0)



include("heat_transfer_plane_geometry.jl")
#include("heat_transfer_spheres.jl")
