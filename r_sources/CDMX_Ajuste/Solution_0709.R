library(pracma)
library(matrixStats)
library(plotly)
library(deSolve)
#
state  <- "CDMX"
ss1    <- 100000
path0  <- "Version_2/"
LL     <- 1000
clevel <- 0.1
#
### Read data ###
pathBase = "~/Insync/sauld@cimat.mx/Google Drive/UNISON/Articles/NovelCovid-19/"
pathSufix = "TwoVaccinesDynamics/multivaccine_multidose_covid_19/r_sources/"
pathData = "CDMX_Ajuste"
path = paste(pathBase, pathSufix, pathData, sep = '')
setwd(path)
Dat1 <- read.csv(paste0(state,"_samples_v16.csv"), head = TRUE)[1:ss1,]
Dat <- Dat1[,-dim(Dat1)[2]]
pIS  <- read.csv(paste0(state, "_qr_1.csv"), head = TRUE)
State_pars <- read.csv("State_parameters.csv", head = TRUE)
### ODE System ###
### Some parameters ###
indx <- which(State_pars$Entities == state)
p_IS <- pchipfun(pIS$Times, as.numeric(pIS$qr))
delta_IS <- 1/State_pars$Nu[indx]
mu_H <- 1/State_pars$Mu[indx]
Dtim <- seq(pIS$Times[1], pIS$Times[length(pIS$Times)], 1)
#
ode1 <- function(t, x1, par1){ 
  S <- x1[1]
  E <- x1[2]
  I_S <- x1[3]
  I_A <- x1[4]
  A <- x1[5]
  H <- x1[6]
  R <- x1[7]
  ##  
  q_H <- par1[1]
  alpha_H <- par1[2]
  delta_IS <- par1[3]
  mu_H <- par1[4]

  alpha_IS <- 1/14    # symptomatic recovery rate
  alpha_IA <- 1/7     # asymptomatic recovery rate
  alpha_A <- 1/10     # ambulatory recovery rate
  q_IS <- 0.85        # % of symptomatic recovered
  q <- 0.45           # asymptomatic infection reduction
  delta_R <- 1/180    # natural immunity
  delta_E <- 1/5.1    # incubation rate
  mu <- 1/(75.5 * 365)  # Natural death rate
  p_E <- 0.2          # % of asymptomatic
  ## Infection force and others
  Ns1 = S + E + I_S + I_A + R + A + H
  Nss1 = Ns1 - (A + H)
  l_f1 = bet(t) * (I_S + q * I_A) / Nss1
  dS = mu * Ns1 - (l_f1 + mu) * S + delta_R * R
  dE = l_f1 * S - (delta_E + mu) * E
  dI_S = (1 - p_E) * delta_E * E -
        (q_IS * alpha_IS + delta_IS * (1 - q_IS) + mu) * I_S
  dI_A = p_E * delta_E * E - (alpha_IA + mu) * I_A
  dA = p_IS(t) * delta_IS * (1 - q_IS) * I_S - (alpha_A + mu)*A
  dH = (1 - p_IS(t)) * delta_IS * (1 - q_IS) * I_S - 
        (q_H * alpha_H + mu_H * (1 - q_H) + mu) * H
  dR = q_IS * alpha_IS * I_S + alpha_IA * I_A + 
        alpha_A * A + q_H * alpha_H * H - (delta_R + mu) * R
  dD = mu_H * (1 - q_H) * H
  dCA = p_IS(t) * delta_IS * (1 - q_IS) * I_S
  dCH = (1 - p_IS(t)) * delta_IS * (1 - q_IS) * I_S
  rhs = list(c(dS, dE, dI_S, dI_A, dA, dH, dR, dD, dCA, dCH))
  return(rhs)
}
### Total population ###
N1 <- State_pars[indx, 9]
### Defining matrix ###
LR <- seq(1, dim(Dat)[1], length.out = LL)
MSol_D <- matrix(NA, nrow = LL, ncol = length(Dtim))
MSol_A <- matrix(NA, nrow = LL, ncol = length(Dtim))
MSol_H <- matrix(NA, nrow = LL, ncol = length(Dtim))
MSol_T <- matrix(NA, nrow = LL, ncol = 10)
#
### Solve system for each sample ###
for (i1 in 1:LL) {
  i2 <- LR[i1]
  bet <- pchipfun(pIS$Times, as.numeric(c(Dat[i2,6], Dat[i2,6:17])))
  q_H <- Dat[i2, 4]
  alpha_H <- 1/Dat[i2, 5]
  E0 <- Dat[i2, 1]
  I_A0 <- Dat[i2, 2]
  I_S0 <- Dat[i2, 3]
  S0 <- N1 - (E0 + I_A0 + I_S0)
  par1 <- c(q_H, alpha_H, delta_IS, mu_H)
  z1 <- c(S0, E0, I_S0, I_A0, 0, 0, 0, 0, 0, 0)
  Y1 <- 
    as.matrix(
      ode(
          func = ode1,
          y = z1,
          times = Dtim,
          parms = par1[1:4],
          method = "rk4")
    )[,1:11]
  MSol_D[i1,] <- c(0, Y1[2:dim(Y1)[1],9] - Y1[1:(dim(Y1)[1] - 1), 9])
  MSol_A[i1,] <- c(0, Y1[2:dim(Y1)[1],10] - Y1[1:(dim(Y1)[1] - 1), 10])
  MSol_H[i1,] <- c(0, Y1[2:dim(Y1)[1],11] - Y1[1:(dim(Y1)[1] - 1), 11])
  MSol_T[i1,] <- Y1[dim(Y1)[1],2:11]
}
#
MF_A <- apply(MSol_A, 2, quantile, probs = c(0.025, 0.5, 0.975))
MF_D <- apply(MSol_D, 2, quantile, probs = c(0.025, 0.5, 0.975))
MF_H <- apply(MSol_H, 2, quantile, probs = c(0.025, 0.5, 0.975))
MF_T <- apply(MSol_T, 2, quantile, probs = c(0.025, 0.5, 0.975))
# Final time quantile solution
### Figures ###
xt <- seq(
    as.Date(pIS$Dates[1],"%d/%m/%y"),
    as.Date(pIS$Dates[length(pIS$Dates)], "%d/%m/%y"),
    1)
dat_T <- read.csv(paste0(state,"_030921.csv"), head = TRUE)
dtemp <- as.Date("2021-03-15") 
# Final plotting date
yind <- which(dat_T$Date == dtemp)
Fig1 <- plot_ly(type = "scatter", mode = "none")
Fig1 <- 
  Fig1 %>% 
    add_trace(
      x = xt,
      y = MF_D[1,],
      mode = "lines",
      line = list(color = 'black', dash = "dash"),
      showlegend = FALSE
    )
Fig1 <- 
  Fig1 %>% 
    add_trace(
      x = xt,
      y = MF_D[2,],
      mode = "lines",
      line = list(color = 'red'),
      showlegend = FALSE
    )
#
Fig1 <- 
  Fig1 %>% 
    add_trace(
      x = xt,
      y = MF_D[3,],
      mode = "lines",
      line = list(color = 'black', dash = "dash"),
      showlegend = FALSE
    )
Fig1 <- 
    Fig1 %>% 
    add_trace(
        x = dat_T$Date[1:yind],
        y = dat_T$Deaths[1:yind],
        type = 'bar', showlegend = FALSE,
        marker = list(
            color = '1E90FF',
            size = 1,
            line = list(color = '1E90FF', width = 1)
            )
    )
Fig1
#
Fig2 <- plot_ly(type = "scatter", mode = "none")
Fig2 <- 
  Fig2 %>% 
    add_trace(
      x = xt,
      y = MF_A[1,],
      mode = "lines",
     line = list(color = 'black', dash = "dash"), showlegend = FALSE
    )
Fig2 <- 
  Fig2 %>% 
    add_trace(
      x = xt,
      y = MF_A[2,], 
      mode = "lines",
      line = list(color = 'red'),
      showlegend = FALSE
    )
Fig2 <- 
  Fig2 %>% 
    add_trace(
      x = xt,
      y = MF_A[3,],
      mode = "lines",
      line = list(
              color = 'black',
              dash = "dash"),
      showlegend = FALSE
    )
Fig2
#
Fig3 <- 
  plot_ly(type = "scatter", mode = "none")
Fig3 <- 
  Fig3 %>% 
    add_trace(
      x = xt,
      y = MF_H[1,],
      mode = "lines",
      line = list(color = 'black', dash = "dash"),
      showlegend = FALSE
    )
Fig3 <- 
    Fig3 %>% 
      add_trace(
        x = xt,
        y = MF_H[2,],
        mode = "lines",
        line = list(color = 'red'),
        showlegend = FALSE
    )
#
Fig3 <- 
  Fig3 %>% 
    add_trace(
      x = xt,
      y = MF_H[3,],
      mode = "lines",
      line = list(color = 'black', dash = "dash"),
      showlegend = FALSE
    )
Fig3
