# This file contains the functions for heat transfer


" Far-field heat transfer between two semi-infinite media,classical approximation "


function farfield_transfer(em1 , em2, T1, T2)

    return sigma*(T2^4-T1^4)/(1.0/em1 + 1.0/em2 - 1.0)

end

function evanescent(ml1 :: MultiLayer, ml2 :: MultiLayer, pol :: Polarization , T1 ,T2)


end
