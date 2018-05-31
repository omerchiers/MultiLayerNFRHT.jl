# This file containes all the optical properties.
# Permitivities are given in the convention : eps_real + i*eps_im

export convert_prop, permittivity, epsi, refr, Drude, Lorentz, Cbn, Sic, Si, Al, Au,Au_latella,Cst

abstract type AbstractMaterial end
abstract type AbstractModel <: AbstractMaterial end



# Model types
struct LorentzModel1{T <: Number} <: AbstractModel
    eps0   :: T
    wp     :: T
    w0     :: T
    gamma0 :: T
end

struct LorentzModel2{T <: Number} <: AbstractModel
    eps0   :: T
    w_lo   :: T
    w_to   :: T
    gamma0 :: T
end

LorentzModel1(eps0,wp,w0,gamma0) = LorentzModel1(eps0,wp,w0,gamma0,0.0)

# Constant permittivity
struct Cst{T <: Number} <: AbstractModel
    val :: T
end
Cst() = Cst(1.0+im*0.0)


"""
    permittivity(material,w)

Compute the dielectric permittivity for a material.

# Arguments
* `material :: AbstractModel` : is a type for the material you want to pass
* `w :: Real` : pulsation in rad s^(-1)

# Example
```julia
julia> permittivity(Sic(),1e13)
6.805820438080644 + 0.002847538251166107im
```
"""
function permittivity()
end


function permittivity(model::LorentzModel1,w) :: Complex128
     eps0   = model.eps0
     wp     = model.wp
     w0     = model.w0
     gamma0 = model.gamma0
     gamma1 = model.gamma1
     return  eps0 + wp^2/(w0*w0 - w*w - im*w*(gamma0+gamma1))
end

Au(gamma1)         = LorentzModel1(9.4,13584.25e12,0.0,109.96e12,gamma1)
Au_latella(gamma1) = LorentzModel1(1.0,1.37e16,0.0,5.32e13,gamma1)
Al                 = LorentzModel1(1.0,2.24e16,0.0,1.22e14)

function permittivity(model::LorentzModel2,w) :: Complex128
     eps0   = model.eps0
     w_lo   = model.wp
     w_to   = model.w0
     gamma0 = model.gamma0
     return  eps0*(w^2-w_lo^2 + im*gamma*w)/(w^2-w_to^2 + im*gamma*w)
end

CbN = LorentzModel2(4.46,2.451e14,1.985e14,9.934e11)
SiC = LorentzModel2(6.7,1.827e14,1.495e14,8.971e11)


permittivity(model::Cst,w) = model.val

Si     = Cst(11.7+im*0.0)
Vacuum = Cst(1.0+im*0.0)

refractiveindex(material :: AbstractMaterial,w) = sqrt(permittivity(material,w))
