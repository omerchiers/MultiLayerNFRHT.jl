# This file containes all the optical properties.
# Permitivities are given in the convention : eps_real - i*eps_im

export convert_prop, permittivity, epsi, refr, Drude, Lorentz, cbn, sic, si, al, au
abstract OptProp

# Generic
immutable Model <: OptProp
    eps0  :: Float64
    wp    :: Float64
    w0    :: Float64
    gamma :: Float64
end

# Dielectrics
immutable Sic <: OptProp end
immutable Cbn <: OptProp end
immutable Si <: OptProp end

# Conductors
immutable Al <: OptProp end
immutable Au <: OptProp end

# Constant permittivity
immutable Cst <: OptProp
    val   :: Complex128
end
Cst() = Cst(1.0+im*0.0)



"""
    permittivity(material,w)

Compute the dielectric permittivity for a material.

# Arguments
* `material :: OptProp` : is a type for the material you want to pass
* `w :: Real` : pulsation in rad s^(-1)

# Example
```julia
julia> permittivity(Sic(),1e13)
6.805820438080644 + 0.002847538251166107im
```
"""
function permittivity() end


function permittivity(material::Model,w) :: Complex128
     eps0  = material.eps0
     wp    = material.wp
     w0    = material.w0
     gamma = material.gamma
    return  eps0 + wp^2/(w0*w0 - w*w - im*w*gamma)
end



function permittivity(material::Cbn,w) :: Complex128
    eps_fin = 4.46 + 0.0*im
    w_lo    = 2.451e14 # rad/s
    w_to    = 1.985e14 # rad/s
    gamma   = 9.934e11 # rad/s

    return eps_fin*(w^2-w_lo^2 + im*gamma*w)/(w^2-w_to^2 + im*gamma*w)
end

function permittivity(material::Sic,w) :: Complex128
    eps_fin = 6.7 + 0.0*im
    w_lo    = 1.827e14 # rad/s
    w_to    = 1.495e14 # rad/s
    gamma   = 8.971e11 # rad/s

    return eps_fin*(w^2-w_lo^2 + im*gamma*w)/(w^2-w_to^2 + im*gamma*w)
end


function permittivity(material::Al,w) :: Complex128
    alum = Model(1.0,2.24e16,0.0,1.22e14)
    return permittivity(alum,w)
end

function permittivity(material::Si,w) :: Complex128
    return 11.7+im*0.0
end

function permittivity(material::Cst,w) :: Complex128
    return material.val
end
