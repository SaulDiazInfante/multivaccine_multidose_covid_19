# Title     : Right hand side definition for the two dosI_S
# vaccination model
# Objective : TODO
# Created by: saul
  # Created on: 9/15/21
# 
source("constructHermiteInterpolationPolynomials.R")
hermitePolynomials <- constructHermiteInterpolationPolynomials()
p_ISV1 <- hermitePolynomials[["p_IS"]]
p_ISV2 <- hermitePolynomials[["p_IS"]]
beta_A <- hermitePolynomials[["beta"]]
p_IS <- hermitePolynomials[["p_IS"]]
#beta_A <- function(t){return(0.05162024)} # For tests
#p_IS <- function(t){return(0.91)}
#p_ISV1 <- function(t){return(0.91)}
#p_ISV2 <- function(t){return(0.91)}
evaluteRhsODEVaccinationModel <- function(t, states, parameters) {
    with(
      as.list(c(states, parameters)), {
        S <- states[1]
        E <- states[2]
        I_S <- states[3]
        I_A <- states[4]
        A <- states[5]
        H <- states[6]
        R <- states[7]
        #
        V_1 <- states[8]
        E_V1 <- states[9]
        I_SV1 <- states[10]
        I_AV1 <- states[11]
        A_V1 <- states[12]
        H_V1 <- states[13]
        #
        V_2 <- states[14]
        E_V2 <- states[15]
        I_SV2 <- states[16]
        I_AV2 <- states[17]
        A_V2 <- states[18]
        H_V2 <- states[19]
        D <- states[20]
        CA <- states[21]
        CH <- states[22]
        X <- states[23]
        #
        Ns <- S + E + I_S + I_A + A + H + R +
          V_1 + V_2 + E_V1 + E_V2 + I_SV1 + I_SV2 +
          I_AV1 + I_AV2 + A_V1 + A_V2 + H_V1 + H_V2
        #
        Nss <- Ns - ( A + H + A_V1 + A_V2 + H_V1 + H_V2)
        lambda_f <- (beta_A(t) / Nss) * 
                (
                  q * (I_A + (1 - p_AV1) * I_AV1 + (1 - p_AV2) * I_AV2) +
                  (I_S + (1 - p_SV1) * I_SV1 + (1 - p_SV2) * I_SV2)
                )
        #
        dS  <- 
          mu * Ns -
          (
            lambda_f + mu + phi_V1 + phi_V2
          ) * S + delta_R * R + gamma_1 * V_1 + gamma_2 * V_2
        dE  <- 
          lambda_f * S  - (delta_E + mu) * E
        dI_S <- 
          (1 - p_E) * delta_E * E - 
          (q_IS * alpha_IS + delta_IS * (1 - q_IS) + mu) * I_S
        dI_A <- 
          p_E * delta_E * E - (alpha_IA + mu) * I_A
        dA <- p_IS(t) * delta_IS * (1 - q_IS) * I_S - (alpha_A + mu) * A
        dH <- 
          (1 - p_IS(t)) * delta_IS * (1 - q_IS) * I_S - 
          (q_H * alpha_H + mu_H * (1 - q_H) + mu) * H
        dR <- 
          q_IS * alpha_IS * I_S + alpha_IA * I_A + alpha_A * A +
          q_H * alpha_H * H + alpha_ISV1 * q_ISV1 * I_SV1 +
          alpha_ISV2 * q_ISV2 * I_SV2 + alpha_IAV1 * I_AV1 +
          alpha_IAV2 * I_AV2 + alpha_AV1 * A_V1 + alpha_AV2 * A_V2 +
          q_HV1 * alpha_HV1 * H_V1 + q_HV2 * alpha_HV2 * H_V2 -
          (delta_R + mu) * R
        dV1 <- phi_V1 * S - ((1 - epsilon_V1) * lambda_f + mu + gamma_1) * V_1
        dEV1 <- (1 - epsilon_V1) * lambda_f * V_1 - (delta_EV1 + mu) * E_V1
        dI_SV1 <- (1 - p_EV1) * delta_EV1 * E_V1 - 
          (q_ISV1 * alpha_ISV1 + delta_ISV1 * (1 - q_ISV1) + mu) * I_SV1
        dI_AV1 <- p_EV1 * delta_EV1 * E_V1 - (alpha_IAV1 + mu) * I_AV1
        dAV1  <- p_ISV1 * delta_ISV1 * (1 - q_ISV1) * I_SV1 -
          (alpha_AV1 + mu) * A_V1
        dHV1 <- (1 - p_ISV1) * delta_ISV1 * (1 - q_ISV1) * I_SV1 -
          (q_HV1 * alpha_HV1 + mu_HV1 * (1 - q_HV1) + mu) * H_V1
        #
        dV2 <- phi_V2 * S - ((1 - epsilon_V2) * lambda_f + mu + gamma_2) * V_2
        dEV2 <- (1 - epsilon_V2) * lambda_f * V_2 - (delta_EV2 + mu) * E_V2
        dI_SV2 <- (1 - p_EV2) * delta_EV2 * E_V2 - 
          (q_ISV2 * alpha_ISV2 + delta_ISV2 * (1 - q_ISV2) + mu) * I_SV2
        dI_AV2 <- p_EV2 * delta_EV2 * E_V2 - (alpha_IAV2 + mu) * I_AV2
        dAV2 <- p_ISV2 * delta_ISV2 * (1 - q_ISV2) * I_SV2 -
          (alpha_AV2 + mu) * A_V2
        dHV2 <- 
          (1 - p_ISV2) * delta_ISV2 * (1 - q_ISV2) * I_SV2 -
          (q_HV2 * alpha_HV2 + mu_HV2 * (1 - q_HV2) + mu) * H_V2
        # 
        dD <- mu_H * (1 - q_H) * H + 
          mu_HV1 * (1 - q_HV1) * H_V1 + mu_HV2 * (1 - q_HV2) * H_V2
        dCA <- p_IS(t) * delta_IS * (1 - q_IS) * I_S
        dCH <- (1 - p_IS(t)) * delta_IS * (1 - q_IS) * I_S
        dX <- (phi_V1 + phi_V2) * (S + E + A + R )
        rhs <- 
          list(
            c(
              dS, dE, dI_S, dI_A, dA, dH, dR, 
              dV1, dEV1, dI_SV1, dI_AV1, dAV1, dHV1,
              dV2, dEV2, dI_SV2, dI_AV2, dAV2, dHV2,
              dD, dCA, dCH, dX 
            )
          )
        return(rhs)
    }
    )
}
