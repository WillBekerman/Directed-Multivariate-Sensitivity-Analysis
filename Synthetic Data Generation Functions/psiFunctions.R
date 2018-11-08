################################################################################
#                         Supporting Psi-functions
################################################################################

psi_huber <- Vectorize( FUN = function(y, Kappa = 2.5) # Huber's psi function (pg 1460 of Rosenbaum 2016)
{
  return(sign(y)*min(c(abs(y), Kappa)))
})

psi_innerTrim <- Vectorize(FUN = function(y, Kappa = 2.5, iota = .5) # Huber's psi function (pg 1460 of Rosenbaum 2016)
{
  assert('0 <= iota < Kappa', {((0 <= iota && iota < Kappa))})
  return(sign(y)*(Kappa / (Kappa - iota)) * max(c(0, min(c(abs(y), Kappa)) - iota)))
})