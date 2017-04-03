# Utility functions
function trapz(x,y)
   n = length(y)
   if n != length(x)
       error("Input x,y must be of same length")
   end
       r = zero(zero(eltype(x))*zero(eltype(y)))
   for i in 2:n
       r += (x[i] - x[i-1]) * (y[i] + y[i-1])
   end
       r/2.0
end
