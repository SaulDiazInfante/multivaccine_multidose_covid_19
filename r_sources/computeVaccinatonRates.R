computeVaccinatonRates <- 
  function(parameters, 
           coverageV1,
           coverageV2,
           administrationTime = 7.0) {
  NSS <- as.numeric(parameters["NSS"])
  phiV1 <- -log(1 - coverageV1 / NSS) / administrationTime
  phiV2 <- -log(1 - coverageV2 / NSS) / administrationTime
  vaccinationRates <- c(phiV1, phiV2)
  return(vaccinationRates)
}