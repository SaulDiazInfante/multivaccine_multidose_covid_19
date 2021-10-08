library(deSolve)
library(rjson)
library(reshape2)
library(dplyr)
library(plotly)
source("loadInitialConditions.R")
source("loadTransferParameters.R")
source("rhsODEVaccinationModel.R")
source("computeVaccinatonRates.R")
#
# two months vaccination. 
# You have to run things so you only stay
# within these values.
#
vaccinated_perweek1 <-
  c(1000, 2000, 1000, 1500, 2000, 7000, 10000, 8000)
vaccinated_perweek2 <-
  c(1000, 2000,1000, 1500, 2000, 7000, 10000, 8000)
Nvacc_schemes <-  length(vaccinated_perweek1)
modelParameters <- fromJSON(file = "ModelParameters.json")
dataFrameModelParameters <- as.data.frame(modelParameters)
timeLine <- seq(0, 55)
parameters_values <- loadTransferParameters(dataFrameModelParameters)
initial_values <- loadInitialConditions(dataFrameModelParameters)
#
splitedTimeLine <-
  split(
    timeLine,
    ceiling(seq_along(timeLine)/7)
  )
# First iteration
currentTimeWeek <- unlist(splitedTimeLine[[1]])
phiV <-
  computeVaccinatonRates(parameters_values, 
                         vaccinated_perweek1[1],
                         vaccinated_perweek2[1])
parameters_values["phi_V1"] <- phiV[1]
parameters_values["phi_V2"] <- phiV[2]
fZero <- evaluteRhsODEVaccinationModel(currentTimeWeek[1], 
                                       initial_values, parameters_values)
currentSolution <- ode(
  y = initial_values,
  times = currentTimeWeek,
  func = evaluteRhsODEVaccinationModel,
  parms = parameters_values
)
currentSolution <- as.data.frame(currentSolution)
odeVarNames <- c("time", "S", "E", "I_S", "I_A", "A", "H", "R",
                 "V_1", "E_V1", "I_SV1", "I_AV1", "A_V1", "H_V1",
                 "V_2",  "E_V2", "I_SV2", "I_AV2", "A_V2", "H_V2",
                 "D", "CA", "CH", "X_k")
colnames(currentSolution) <- odeVarNames
#
for (k in 2:length(splitedTimeLine)) {
  # Select interval of integration
  currentTimeWeek <- unlist(splitedTimeLine[[k]])
  # Update parameters
  initial_values <- tail(currentSolution[2:24], n = 1)
  phiV <- 
    computeVaccinatonRates(parameters_values, 
                           vaccinated_perweek1[k],
                           vaccinated_perweek2[k])
  parameters_values["phi_V1"] <- phiV[1]
  parameters_values["phi_V2"] <- phiV[2]
  newSolution <- ode(
    y = as.numeric(initial_values),
    times = currentTimeWeek,
    func = evaluteRhsODEVaccinationModel,
    parms = parameters_values
  )
  newSolution <- as.data.frame(newSolution)
  colnames(newSolution) <- odeVarNames
  # Joint solution
  currentSolution <- bind_rows(currentSolution, newSolution)
}
#
#
# required by ggplot: data object must be a data frame
initial_values <- loadInitialConditions(dataFrameModelParameters)
N_0 <- sum(initial_values)- (421363.8 + 49450.24)
N_t <-
  (
    currentSolution$S + currentSolution$E + currentSolution$I_S + 
      currentSolution$I_A + currentSolution$A + currentSolution$H + 
      currentSolution$R + currentSolution$V_1 + currentSolution$E_V1 + 
      currentSolution$I_SV1 + currentSolution$I_AV1 + currentSolution$A_V1 + 
      currentSolution$H_V1 + currentSolution$V_2 + currentSolution$E_V2 + 
      currentSolution$I_SV2 + currentSolution$I_AV2 + currentSolution$A_V2 + 
      currentSolution$H_V2 + currentSolution$D
  ) / N_0

fig00 <- ggplot(currentSolution, aes(x = time, y = N_t))
fig00 <-  fig00 + geom_line()
plotly::ggplotly(fig00)

figVac <- ggplot(currentSolution, aes(x = time, y = X_k))
figVac <- figVac + geom_line()
plotly::ggplotly(figVac)

figS <- plot_ly(currentSolution, x = ~time, 
                 y = ~S, 
                 type = 'bar')

figI_S <- plot_ly(currentSolution, x = ~time, 
                y = ~I_S, 
                mode = 'lines+markers')

figI_A <- plot_ly(currentSolution, x = ~time, 
                  y = ~I_A, 
                  mode = 'lines+markers')
