

" Transfer function between two spheres with identical radius and radius << λ "
function transfer_spheres(d,radius,eps1 ,eps2,w)

    trans1 = 9.0*c0^2/(2.0*w^2*d^2) + 9.0*c0^4/(2.0*w^4*d^4) + 27*c0^6/(2.0*w^6*d^6)
    trans2 = 9.0*c0^2/(2.0*w^2*d^2) + 9.0*c0^4/(2.0*w^4*d^4)

    trans = 0.0+0.0*im

    for pol1 in (te(),tm())
        for pol2 in (te(),tm())
            if pol1==pol2
                tran = (real(t1(pol1,radius,eps1,w)) + abs(t1(pol1,radius,eps1,w))^2)*(real(t1(pol1,radius,eps2,w)) + abs(t1(pol1,radius,eps2,w))^2)*trans1
            else
                tran = (real(t1(pol1,radius,eps1,w)) + abs(t1(pol1,radius,eps1,w))^2)*(real(t1(pol2,radius,eps2,w)) + abs(t1(pol2,radius,eps2,w))^2)*trans2
            end
            trans += tran
        end
    end
    return (2.0/pi)*trans
end


" Mie coefficient a1 "
function t1(pol :: te(),radius,eps,w ; mu ::Float64 = 1.0)
    rst = radius*w/c0
    a1 = im*2.0*(rst^3)*(eps-1.0)/3.0/(eps+2.0)
    a2 = 2.0*im*(2.0-3.0*eps+eps^2*(1.0+mu))/5.0/(2.0+eps)^2
    a3 = -4.0*(rst^6)*(eps-1.0)^2/9.0/(2.0+eps)
    return a1+a2+a3
end

" Mie coefficient b1 "
function t1(pol :: tm(),radius,eps,w ; mu ::Float64 = 1.0)
    rst = radius*w/c0
    a1 = im*2.0*(rst^3)*(mu-1.0)/3.0/(mu+2.0)
    a2 = 2.0*im*(2.0-3.0*mu+mu^2*(1.0+eps))/5.0/(2.0+mu)^2
    a3 = -4.0*(rst^6)*(mu-1.0)^2/9.0/(2.0+mu)
    return a1+a2+a3
end


" Spectral heat transfer between two identical spheres with radius << λ "

# function heat_transfer_spheres_w(d,radius,eps1,eps2,w)
#
#
# end
