## Aux parameters ##
library(pracma)
library(tidyverse)
library(dplyr)
constructHermiteInterpolationPolynomials <- 
  function(fileName = "cdmxInterpolationData.csv"){
  pIS <- as.data.frame(read.csv(fileName, head = TRUE))
  p_IS <- 
    pchipfun(
      pIS$Times[1:length(pIS$Times)],
      as.numeric(pIS$qr[1:length(pIS$Times)])
    )
  beta <- 
    pchipfun(
      pIS$Times[1:length(pIS$Times)],
      as.numeric(pIS$beta_ti[1:length(pIS$beta_ti)])
    )
  p_ISV1 <- pchipfun(
    pIS$Times[1:length(pIS$Times)],
    as.numeric(pIS$qr[1:length(pIS$Times)])
  )
    
  p_ISV2 <-     pchipfun(
    pIS$Times[1:length(pIS$Times)],
    as.numeric(pIS$qr[1:length(pIS$Times)])
  )
  hermitePolynomials <-
    list(p_IS = p_IS, beta = beta, p_ISV1 = p_ISV1, p_ISV2 = p_ISV2)
  #names(hermitePolynomials) <- c("p_IS", "beta", "p_ISV1", "p_ISV2")
  return(hermitePolynomials)
} 
