module Tmp
 include("optical_properties.jl")
 include("fresnel_coefficients.jl")
end


function abeles_matrix(pol::te,layer :: Layer )
    return m
end

function abeles_matrix(pol::tm,layer :: Layer )
    return m
end


function rt(struct :: Ml , pol::te, kz,  mlayer :: Vector{Layer})
    ## Implement Abeles formalism
    m =ones(2)
    for i=1:length(mlayer)
        m = m*abeles_matrix(pol,mlayer[i])
    end
    r = m[1,1]
    t = m[2,2]
end

function rt(struct :: Ml, pol::tm, kz,  mlayer :: Vector{Layer})
    ## Implement Abeles formalism
    m =ones(2)
    for i=1:length(mlayer)
        m = m*abeles_matrix(pol,mlayer[i])
    end
    r = m[1,1]
    t = m[2,2]
end


function convert_prop(::refr,refrind)
# Receives the refractive index and returns the permittivity
    epsr  =
    epsim =
    return epsr + im*epsim
end

function permittivity(material::Drude,w; wp= , gamma =   ,eps0 =)
    eps0 - wp.*wp./(w.*(w. + im*gamma))
end

function permittivity(material::Lorentz,w; w0= , gamma =   ,eps0 =)
    eps0 - wp.*wp./((w.-w0).^2 + im.*w.*gamma))
end


function convert_prop(eps::epsi,w)
# Receives the permittivity and returns the refractive index

    epsr = real(eps)
    epsi = imag(eps)
    n = sqrt(0.5*(epsr+sqrt(epsr^2.0+epsi^2.0)))
    k = sqrt(0.5*(-epsr+sqrt(epsr^2.0+epsi^2.0)))
    return n+im*k
end

## New type organization for Bulk
abstract Structure
abstract OptProp

immutable Layer <: Structure
    material  :: OptProp
    thickness :: Float64 # in meter. For semi-infinite media thickness = 0
end

Layer(material) = Layer(material,0.0) # Defines a semi-infinite layer

typealias MultiLayer Vector{Layer}  # Defines a multilayer
