rm(list=ls())
library(deSolve)
basic_model <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
    Ns <- S + E + IS + IA + R + A + H + L + V1 + V2 + EV1 + EV2 +ISV1 + IAV1 
    Nss <- Ns - A - H - AV1 - AV2 - HV1 - HV2
    l_f <- 
      (beta/Nss)*(q*(IA+(1-pAV1)*IAV1+(1-pAV2)*IAV2)+ (IS+(1-pSV1)*ISV1+(1-pSV2)*ISV2))
    
    dS  <- mu*Nb - (l_f + mu + fi_V1 + fi_V2)*S + d_R*R + ga_1*V1 + ga_2*V2 
    dE  <- l_f *S  - (d_E *E + mu)*E
    dIS <- (1-p_E)*d_E *E - (q_IS*a_IS +d_IS*(1-q_IS) + mu)* IS
    dIA <- p_E*d_E*E - (a_IA+mu)*IA
    dA  <- p_IS*d_IS*(1-q_IS)*IS-(a_A+mu)*A
    dH  <- (1-p_IS)*d_IS*(1-q_IS)*I_S-(q_H*a_H+mu_H*(1-q_H)+mu)*H
    dR  <- q_IS*a_IS*IS + a_IA*IA + a_A*A + q_H*a_H*H + a_ISV1*q_ISV1*ISV1 + a_ISV2*q_ISV2*ISV2 + a_IAV1*IAV1 + a_IAV2*IAV2 + a_AV1*AV1 + a_AV2*AV2 + q_HV1*a_HV1*HV1 + q_HV2*a_HV2*HV2 - (d_R + mu)*R
    dL  <- d_S*S-(1-e_L)*l_f*L-(d_L+mu)*L-(fi_V1+fi_V2)*L
      
    dV1   <- fi_V1*S - ((1-e_V1)*l_f + mu +ga_1)*V1
    dEV1  <- (1-e_V1)*l_f*V1 - (d_EV1 + mu)*EV1
    dISV1 <- (1-p_EV1)*d_EV1 *EV1 - (q_ISV1*a_ISV1 +d_ISV1*(1-q_ISV1) + mu)* ISV1
    dIAV1 <- p_EV1*d_EV1*EV1 - (a_IAV1+mu)*IAV1
    dAV1  <- p_ISV1*d_ISV1*(1-q_ISV1)*ISV1-(a_AV1+mu)*AV1
    dHV1  <- (1-p_ISV1)*d_ISV1*(1-q_ISV1)*I_SV1-(q_HV1*a_HV1+mu_HV1*(1-q_HV1)+mu)*HV1
    
    dV2   <- fi_V2*S - ((1-e_V2)*l_f + mu +ga_1)*V2   
    dEV2  <- (1-e_V2)*l_f*V2 - (d_EV2 + mu)*EV2
    dISV2 <- (1-p_EV2)*d_EV2 *EV2 - (q_ISV2*a_ISV2 +d_ISV2*(1-q_ISV2) + mu)* ISV2
    dIAV2 <- p_EV2*d_EV2*EV2 - (a_IAV2+mu)*IAV2
    dAV2  <- p_ISV2*d_ISV2*(1-q_ISV2)*ISV2-(a_AV2+mu)*AV2
    dHV2  <- (1-p_ISV2)*d_ISV2*(1-q_ISV2)*I_SV2-(q_HV2*a_HV2+mu_HV2*(1-q_HV2)+mu)*HV2   
    
    dD  <- mu_H*(1-q_H)*H + mu_HV1*(1-q_HHV1)*HV1 + mu_HV2*(1-q_HHV2)*HV2
    dCA <- p_IS*d_IS*(1-q_IS)*IS
    dCH <- (1-p_IS)*d_IS*(1-q_IS)*IS
    
    return(list(c(dS, dE, dIS, dIA, dA, dH, dR, dL , dV1 , dEV1, dISV1 , dIAV1 , dHV1 , dV2 , dEV2, dISV2 , dIAV2 , dHV2 , dD, dCA, dCH)))
  })
}
#
N <- 10000
a_A = 1
a_H <-1
parameters_values <- 
  c(
    a_IS=1/14,
    a_IA=1/7,
    a_A=1/10, 
    a_AV1=a_A,
    a_AV2=a_A,
    a_H=1, 
    a_HV1=a_H,
    a_HV2=a_H,
    a_IAV1=1,
    a_IAV2=1,
    a_ISV1=1,
    a_ISV2=1,
    b_A=1,
    b_S=1,
    d_R=1/180,
    d_E=1/5.1,
    d_IS=1,
    d_S=1,
    d_ISV1=1,
    d_ISV2=1,
    e_L=1,
    e_V1=1,
    e_V2=1,
    fi_V1=1,
    fi_V2=1,
    ga_1=1,
    ga_2=1,
    mu=1/(75.5*365),
    mu_HV1=1,
    mu_HV2=1,
    p_E=0.2, 
    p_EV1=1,
    p_EV2=1,
    p_IS=1,
    p_ISV1=1,
    p_ISV2=1,
    q_IS=0.85,
    q=1,
    q_ISV1=1,  
    q_ISV2=1,
    q_HV1=1,
    q_HV2=1
)
initial_values <- c(
  S=N-1, 
  E=0, 
  IS=0, 
  IA=0, 
  A=0, 
  H=0, 
  R=0, 
  D=0,
  L=0
)
# Indicamos el n? de d?as a simular
time_values <- seq(0, 60)
# Simulamos
sir_values_1 <- ode(
  y = initial_values,
  times = time_values,
  func = basic_model,
  parms = parameters_values
)


sir_values_1 <- as.data.frame(sir_values_1)
with(sir_values_1, {
  plot(time, I, type = "l", col = "blue",
       xlab = "tiempo (dias)", ylab = "numero de personas")
  lines(time, S, col = "red")
  lines(time, R, col = "green")
  #xlim=c(0, 60), ylim=c(0, 1000))
})
legend("right", c("infecciosos", "susceptibles", "recuperados"),
       col = c("blue", "red", "green"), lty = 1, bty = "n")
