computeVaccinatonRates <- 
  function(parameters, 
           coverageV1,
           coverageV2,
           administrationTime = 7.0) {
  NSS <- as.numeric(parameters["NSS"])
  phiV1 <- -log(NSS - coverageV1) / administrationTime
  phiV2 <- -log(NSS - coverageV2) / administrationTime
  vaccinationRates <- c(phiV1, phiV2)
  return(vaccinationRates)
}