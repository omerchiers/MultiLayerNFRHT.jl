# This file defines the reflection and transmission coefficients of a multilayer

using StaticArrays

abstract Polarization
immutable te <: Polarization end
immutable tm <: Polarization end


abstract Structure

immutable Bulk{T <:OptProp,U <:OptProp} <: Structure
    ep1 :: T
    ep2 :: U
end
# Outer constructor for the Bulk type when one medium is vacuum
Bulk(ep2) = Bulk(Cst(1.0+im*0.0),ep2)

immutable Layer <: Structure
    material  :: OptProp
    thickness :: Float64 # in meter. For semi-infinite media thickness = 0
end
Layer(material)  = Layer(material,0.0)


#typealias MultiLayer Union{Vector{Layer},Vector{Layer{Cst}},Vector{Layer{Al}}}
typealias MultiLayer Vector{Layer}



function compute_kz(kx,eps,w)
    eps = conj(eps)
    sqrt(eps*(w/c0)^2 - kx*kx)
end

function compute_kx(theta_0,eps,w)
    eps=conj(eps)
    sqrt(eps)*sin(theta_0)*w/c0
end

function compute_ang(kx,eps,w)
    asin(real(kx)/real(sqrt(eps))/(w/c0))
end



# Fresnel coefficient of a semi-infinite medium : perpendicular
function rt(pol :: te, eps1,eps2, k0z,k2z,w)
  # k0z = incident wavevector from medium 1
    eps1 = conj(eps1)
    eps2 = conj(eps2)

    r = (k0z - k2z)/(k0z + k2z)
    t = (2.0+im*0.0)*k0z/(k0z + k2z)
    return r,t
end

# Fresnel coefficient of a semi-infinite medium : parallel
function rt(pol :: tm, eps1,eps2, k0z,k2z,w)
  # k0z = incident wavevector from medium 1
    eps1 = conj(eps1)
    eps2 = conj(eps2)

    r = (eps2*k0z-eps1*k2z)/(eps2*k0z+eps1*k2z)
    t = (2.0+im*0.0)*sqrt(eps2)*sqrt(eps1)*k0z/(eps2*k0z+eps1*k2z)
    return r,t
end

# Fresnel coefficient of a semi-infinite medium : as a function of incidence angle and frequency
function rt(struct :: Bulk, pol :: Polarization, kx ,w)

    eps1 = permittivity(struct.ep1,w)
    eps2 = permittivity(struct.ep2,w)

    k0z   = compute_kz(kx,eps1,w)
    k2z   = compute_kz(kx,eps2,w)

    r,t=rt(pol,eps1,eps2,k0z,k2z,w)
    return r, t
end


# Fresnel coefficient of a multilayered semi-infinite medium : as a function of incidence angle and frequency

function rt(struct :: MultiLayer, pol :: Polarization, kx ,w)
    ab_matrix  =[1.0+0.0*im  0.0+0.0*im ;
              0.0+0.0*im  1.0+0.0*im ]
    if real(kx)<w/c0 || real(kx)>w/c0
        ab_matrix = abeles_matrix!(ab_matrix,struct,pol,kx,w)
        return ab_matrix[2,1]/ab_matrix[1,1], (1.0+im*0.0)/ab_matrix[2,1]
    else
        return 1.0+0.0*im, 0.0+0.0*im
    end
end

function abeles_matrix!(ab_matrix,struct :: MultiLayer, pol :: Polarization, kx, w)

        for i=2:length(struct)
            eps1 :: Complex128 = permittivity(struct[i-1].material,w)
            eps2 :: Complex128 = permittivity(struct[i].material,w)
            k0z   = compute_kz(kx,eps1,w)
            k2z   = compute_kz(kx,eps2,w)
            r,t = rt(pol, eps1,eps2, k0z,k2z,w)

            interface_matrix  = [1.0+im*0.0  r ; r  1.0+im*0.0 ]./t
            layer_matrix      = [exp(im*struct[i].thickness*k2z)  0.0+0.0*im ;
                                0.0+0.0*im  exp(-im*struct[i].thickness*k2z) ]
            ab_matrix  = ab_matrix*(interface_matrix*layer_matrix)
        end
    return ab_matrix

end
