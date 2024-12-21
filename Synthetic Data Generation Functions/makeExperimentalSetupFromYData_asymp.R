makeExperimentalSetupFromYData_asymp <-function(YData, psi)
{
  omega = YData$omega_k
  YOverOmega = sweep(YData$YData, 2, omega, "/")
  psi_YOverOmega = apply(YOverOmega, MARGIN = 2, FUN = psi)
  return(list(psi_YOverOmega = psi_YOverOmega))
}