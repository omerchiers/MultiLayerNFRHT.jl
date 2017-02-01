
using Plots

include("physical_constants.jl")
include("optical_properties.jl")
include("rt_coefficients.jl")
include("emissivity.jl")
include("heat_transfer.jl")


function main(material ::BulkOrMultiLayer ,w1,w2,N)

    wv    = collect(linspace(w1,w2,N))
    em_w  = zeros(Float64,N)
#    eps_w = zeros(Complex128,N)

    for i=1:N
        em_w[i] = emissivity_w( material, wv[i])
#        eps_w[i] = permittivity(sic(),wv[i])
    end
    println(mean(em_w))
    #plot(wv,real(eps_w))
    plot(wv,em_w)
    xaxis!(:log10)
end
