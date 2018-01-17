
abstract type TotalField end
struct Evanescent  <: TotalField end
struct Propagative <: TotalField end


include("heat_transfer_plane_geometry.jl")
