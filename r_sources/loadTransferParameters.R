# Title     : TODO
# Objective : TODO
# Created by: saul
# Created on: 9/15/21
library(dplyr)
loadTransferParameters <- function(parametersDataFrame){
  transferParametersNames <-  dataFrameModelParameters %>%
    select(!starts_with("initialConditions")) %>% names()
  transferParameters <- dataFrameModelParameters %>%
    select(.dots = transferParametersNames)
  colnames(transferParameters) <- transferParametersNames
  return(transferParameters)
}
