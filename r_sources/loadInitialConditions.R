# Title     : TODO
# Objective : TODO
# Created by: saul
# Created on: 9/15/21
# Title     : TODO
# Objective : TODO
# Created by: saul
# Created on: 9/15/21
library(dplyr)
loadInitialConditions <- function(parametersDataFrame){
  initialConditionNames <-  parametersDataFrame %>%
    select(ends_with("0")) %>% names()
  initialCondition <- parametersDataFrame %>%
    select(.dots=initialConditionNames)
  initialCondition_var <- as.vector(t(initialCondition))
  odeVarNames <- c("S", "E", "I_S", "I_A", "A", "H", "R",
                   "V_1", "E_V1", "I_SV1", "I_AV1", "A_V1", "H_V1",
                   "V_2",  "E_V2", "I_SV2", "I_AV2", "A_V2", "H_V2",
                   "D", "CA", "CH", "X_k")
    colnames(initialCondition) <- odeVarNames
    ans <- initialCondition
  return (as.vector(t(ans)))
}

