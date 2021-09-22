#In this model we start vaccination on Feb 14, 2021.
#Uses some outputs from Solution_0709.R by A. Acuña

#This code includes vaccination. Vaccination has a weekly strategy.
#Each week, a number N_k individuals from the S class get vaccinated.
#The scale of the problem is in days.

#Information obtained from Acuña Code

#Information obtained from (paper)

rm(list=ls())
# Cargamos paquete deSolve
library(deSolve)

#two months vaccination. You have to run things so you only stay 
#within these values.
#IMPORTANT: Simulation for 60 days at most.
#vaccinated_perweek1<-c(100000,200000,100000,150000,20000,70000,10000,8000,4000,400) #TBD
#vaccinated_perweek2<-c(100000,200000,100000,150000,20000,70000,10000,8000,4000,400) #TBD

vaccinated_perweek1<-replicate(53, 0) #TBD
vaccinated_perweek2<-replicate(53, 0) #TBD


Nvacc_schemes<-length(vaccinated_perweek1)


# Indicamos par?metros iniciales

parameters_values <- c(   
  alpha_IS=1/14, # symptomatic recovery rate
  alpha_IA=1/7, # asymptomatic recovery rate
  alpha_A=0.1,#1/10, # ambulatory recovery rate
  alpha_AV1=0.1,#alpha_A,
  alpha_AV2=0.1,#alpha_A,
  alpha_H=0.1054311, #1/9.484865, #Solution_0709.R 
  alpha_HV1=0.1054311,#alpha_H,
  alpha_HV2=0.1054311,#alpha_H,
  alpha_IAV1=1/7,
  alpha_IAV2=1/7,
  alpha_ISV1=1/14,
  alpha_ISV2=1/14,
  beta=0.05162024, #Solution_0709.R
  delta_R=1/180, # natural immunity
  delta_E=1/5.1, # incubation rate 
  delta_EV1=1/5.1,
  delta_EV2=1/5.1,
  delta_IS=1/3.5, #0.285714285714286
  delta_ISV1=1/3.5,
  delta_ISV2=1/3.5,
  epsilon_V1=0.05,   #pfizer
  epsilon_V2=0.24,   #AstraZeneca 
  gamma_1=1/180,
  gamma_2=1/180,
  mu=1/(75.5*365), #listo # Natural death rate
  mu_H=1/7, #0.142857142857143
  mu_HV1=1/7,
  mu_HV2=1/7,
  p_AV1=1,
  p_AV2=1,
  p_E=0.2, # % of asymptomatic
  p_EV1=0.2,
  p_EV2=0.2,
  p_IS=0.91, #este depende de t toma el ultimo valor como beta
  p_ISV1=0.91, #igual
  p_ISV2=0.91, #igual
  p_SV1=1,
  p_SV2=1,
  q_IS=0.85, # % of symptomatic recovered
  q=0.45,
  q_ISV1=0.85,  
  q_ISV2=0.85,
  q_H= 0.6194573, #Solution_0709.R
  q_HV1=0.6194573,
  q_HV2=0.6194573
)
initial_values <- c(
  S=8034203,  
  E=13964.548, 
  IS=39595.74, 
  IA=6120.975,
  A=23866.84, 
  H=1656.947, 
  R=877310.1, 
  V1=0,
  EV1=0,
  ISV1=0,
  IAV1=0,
  AV1=0,
  HV1=0,
  V2=0,
  EV2=0,
  ISV2=0,
  IAV2=0,
  AV2=0,
  HV2=0,
  D=21479.00,
  CA=421363.8,
  CH=49450.24
)
N_0<-sum(initial_values)-421363.8-49450.24
#N<- S+E+IS+IA+A+H+R+V1+EV1+ISV1+IAV1+AV1+HV1+V2+EV2+ISV2+IAV2+AV2+HV2+D
# Indicamos el n? de d?as a simular
time_values <- seq(0, 365)
actual_vacc_scheme=1 #indicates the vaccination week.
S_k=8034203 #This is the population at the beginning of the kth week

# Indicamos las ecuaciones diferenciales del modelo SIR
basic_model <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {

    scheme_num=time/7
    if (scheme_num < actual_vacc_scheme) {
      phi_V1<- -(1/7)*log((S_k-vaccinated_perweek1[actual_vacc_scheme])/(S_k))
      phi_V2<- -(1/7)*log((S_k-vaccinated_perweek2[actual_vacc_scheme])/(S_k))
    }else{
      actual_vacc_scheme<<-actual_vacc_scheme+1 
      S_k<<-S 
      phi_V1=-(1/7)*log((S_k-vaccinated_perweek1[actual_vacc_scheme])/(S_k))
      phi_V2=-(1/7)*log((S_k-vaccinated_perweek2[actual_vacc_scheme])/(S_k))
    }
    
    Ns <- S + E + IS + IA + A + H + R + V1 + V2 + EV1 + EV2 +ISV1 + ISV2 + IAV1 +IAV2 + AV1 + AV2 + HV1 + HV2 
    Nss <- Ns - A - H - AV1 - AV2 - HV1 - HV2
    l_f <- (beta/Nss)*(q*(IA+(1-p_AV1)*IAV1+(1-p_AV2)*IAV2)+ (IS+(1-p_SV1)*ISV1+(1-p_SV2)*ISV2))
    
    dS  <- mu*Ns - (l_f + mu + phi_V1 + phi_V2)*S + delta_R*R + gamma_1*V1 + gamma_2*V2 
    dE  <- l_f *S  - (delta_E *E + mu)*E
    dIS <- (1-p_E)*delta_E *E - (q_IS*alpha_IS +delta_IS*(1-q_IS) + mu)* IS
    dIA <- p_E*delta_E*E - (alpha_IA+mu)*IA
    dA  <- p_IS*delta_IS*(1-q_IS)*IS-(alpha_A+mu)*A
    dH  <- (1-p_IS)*delta_IS*(1-q_IS)*IS-(q_H*alpha_H+mu_H*(1-q_H)+mu)*H
    dR  <- q_IS*alpha_IS*IS + alpha_IA*IA + alpha_A*A + q_H*alpha_H*H + alpha_ISV1*q_ISV1*ISV1 + alpha_ISV2*q_ISV2*ISV2 + alpha_IAV1*IAV1 + alpha_IAV2*IAV2 + alpha_AV1*AV1 + alpha_AV2*AV2 + q_HV1*alpha_HV1*HV1 + q_HV2*alpha_HV2*HV2 - (delta_R + mu)*R
    
    dV1   <- phi_V1*S - ((1-epsilon_V1)*l_f + mu +gamma_1)*V1
    dEV1  <- (1-epsilon_V1)*l_f*V1 - (delta_EV1 + mu)*EV1
    dISV1 <- (1-p_EV1)*delta_EV1 *EV1 - (q_ISV1*alpha_ISV1 +delta_ISV1*(1-q_ISV1) + mu)* ISV1
    dIAV1 <- p_EV1*delta_EV1*EV1 - (alpha_IAV1+mu)*IAV1
    dAV1  <- p_ISV1*delta_ISV1*(1-q_ISV1)*ISV1-(alpha_AV1+mu)*AV1
    dHV1  <- (1-p_ISV1)*delta_ISV1*(1-q_ISV1)*ISV1-(q_HV1*alpha_HV1+mu_HV1*(1-q_HV1)+mu)*HV1
    
    dV2   <- phi_V2*S - ((1-epsilon_V2)*l_f + mu +gamma_2)*V2   
    dEV2  <- (1-epsilon_V2)*l_f*V2 - (delta_EV2 + mu)*EV2
    dISV2 <- (1-p_EV2)*delta_EV2 *EV2 - (q_ISV2*alpha_ISV2 +delta_ISV2*(1-q_ISV2) + mu)* ISV2
    dIAV2 <- p_EV2*delta_EV2*EV2 - (alpha_IAV2+mu)*IAV2
    dAV2  <- p_ISV2*delta_ISV2*(1-q_ISV2)*ISV2-(alpha_AV2+mu)*AV2
    dHV2  <- (1-p_ISV2)*delta_ISV2*(1-q_ISV2)*ISV2-(q_HV2*alpha_HV2+mu_HV2*(1-q_HV2)+mu)*HV2   
    
    dD  <- mu_H*(1-q_H)*H + mu_HV1*(1-q_HV1)*HV1 + mu_HV2*(1-q_HV2)*HV2
    dCA <- p_IS*delta_IS*(1-q_IS)*IS
    dCH <- (1-p_IS)*delta_IS*(1-q_IS)*IS
    
    return(list(c(dS,dE,dIS,dIA,dA,dH,dR,dV1,dEV1,dISV1,dIAV1,dAV1,dHV1,dV2,dEV2,dISV2,dIAV2,dAV2,dHV2,dD,dCA,dCH)))
  })
}

# Simulamos
salida <- ode(
  y = initial_values,
  times = time_values,
  func = basic_model,
  parms = parameters_values
)


## Ploting
salida.df = as.data.frame(salida) # required by ggplot: data object must be a data frame
library(reshape2)

N_t<-(salida.df$S+salida.df$E+salida.df$IS+salida.df$IA+salida.df$A+salida.df$H+salida.df$R+salida.df$V1+salida.df$EV1+salida.df$ISV1+salida.df$IAV1+salida.df$AV1+salida.df$HV1+salida.df$V2+salida.df$EV2+salida.df$ISV2+salida.df$IAV2+salida.df$AV2+salida.df$HV2+salida.df$D)/N_0
plot(N_t)
#


# Creamos el gr?fico 1
#salida <- as.data.frame(salida)
#with(salida, {
#  plot(time, IA, type = "l", col = "blue", xlab = "tiempo (dias)", ylab = "numero de personas")
#  # lines(time, S, col = "red")
#  # lines(time, D, col = "green")
#})
#legend("right", c("infecciosos", "susceptibles", "recuperados"),
#       col = c("blue", "red", "green"), lty = 1, bty = "n")


#salida.m = melt(salida.df, id.vars='time') # this makes plotting easier by puting all variables in a single column
#library(ggplot2)
#p <- ggplot(salida.m, aes(time, value, color = variable)) + geom_point()
#print(p)   



