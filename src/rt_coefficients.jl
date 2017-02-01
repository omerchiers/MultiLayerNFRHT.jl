# This file defines the reflection and transmission coefficients of a multilayer


abstract Polarization
immutable te <: Polarization end
immutable tm <: Polarization end


abstract Structure

immutable Bulk <: Structure
    ep1 :: OptProp
    ep2 :: OptProp
end
# Outer constructor for the Bulk type when one medium is vacuum
Bulk(ep) = Bulk(Cst(1.0+im*0.0),ep)

immutable Layer <: Structure
    material  :: OptProp
    thickness :: Float64 # in meter. For semi-infinite media thickness = 0
end
Layer(material) = Layer(material,0.0)

typealias MultiLayer Vector{Layer}


function compute_kz(kx,eps,w)
    sqrt(eps*(w/c0)^2 - kx*kx)
end

function compute_kx(theta_0,eps,w)
    sqrt(eps)*sin(theta_0)*w/c0
end

function compute_ang(kx,eps,w)
    asin(real(kx)/real(sqrt(eps))/(w/c0))
end



# Fresnel coefficient of a semi-infinite medium : perpendicular
function rt(struct :: Bulk, pol :: te, k0z,k2z,w)
  # k0z = incident wavevector from medium 1
    t = (2.0+im*0.0)*k0z/(k0z + k2z)
    r = (k0z - k2z)/(k0z + k2z)
    return r,t
end

# Fresnel coefficient of a semi-infinite medium : parallel
function rt(struct :: Bulk, pol :: tm, k0z,k2z,w)
  # k0z = incident wavevector from medium 1

    eps1 =permittivity(struct.ep1,w)
    eps2 =permittivity(struct.ep2,w)

    t = (2.0+im*0.0)*eps2*k0z/(eps2*k0z+eps1*k2z)
    r = (eps2*k0z-eps1*k2z)/(eps2*k0z+eps1*k2z)
    return r,t
end

# Fresnel coefficient of a semi-infinite medium : as a function of incidence angle and frequency
function rt(struct :: Bulk, pol :: Polarization, kx ,w)

    eps1 = permittivity(struct.ep1,w)
    eps2 = permittivity(struct.ep2,w)

    k0z   = compute_kz(kx,eps1,w)
    k2z   = compute_kz(kx,eps2,w)

    r,t=rt(struct, pol, k0z,k2z,w)
    return r, t
end


# Fresnel coefficient of a multilayered semi-infinite medium : as a function of incidence angle and frequency

function rt(struct :: MultiLayer, pol :: Polarization, kx ,w)

    ab_matrix  = abeles_matrix(struct,pol,kx,w)

    r = ab_matrix[2,1]/ab_matrix[1,1]
    t = (1.0+im*0.0)/ab_matrix[1,1]

    return r, t
end

function abeles_matrix(struct :: MultiLayer, pol :: Polarization, kx, w) :: Array{Complex128,2}
    ab_matrix  = [1.0+0.0*im  0.0+0.0*im ; 0.0+0.0*im  1.0+0.0*im ]

    for i=2:length(struct)
        mat1  = struct[i-1].material
        mat2  = struct[i].material

        eps2  = permittivity(mat2,w)
        k2z   = compute_kz(kx,eps2,w)

        r,t = rt(Bulk(mat1,mat2), pol , kx ,w)

        interface_matrix :: Array{Complex128,2} = [1.0+im*0.0  r ; r  1.0+im*0.0 ]/t
        beta             = struct[i].thickness*k2z
        layer_matrix     :: Array{Complex128,2} = [exp(-im*beta) 0.0+im*0.0 ; 0.0+im*0.0 exp(im*beta)]
        ab_matrix        = ab_matrix*(interface_matrix*layer_matrix)
    end
    return ab_matrix

end
