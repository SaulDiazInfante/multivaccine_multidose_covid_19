import numpy as np
from pytwalk3 import pytwalk3 as pytwalk
import pandas as pd
from plots_setup import *
from scipy.interpolate import pchip
from analysis_v10 import getMAP, getsample

# Plotting parameters
set_plot_parameters()

entity = "CDMX"
period = "daily"
#
# MCMC run
#
version = "v16"
data_path    = ""
samples_path = ""
figures_path = ""
ni = 10000  # 4000000
burnin = 1000  # 100000
ess = 1000  # 100000  # effective sample size
samples_filename = samples_path + entity + '_samples_' + version + '.csv'
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Loading observed data
data_filename = data_path + entity + "_covid_" + period + "_movav.csv"
dates_filename = samples_path + entity + "_keydates.csv"
pop_filename = data_path + "codigos_estados.csv"
qr_filename = data_path + entity + "_qr_1.csv"
Raw = pd.read_csv(data_filename, header=0, sep=",")  # raw data
Keydates = pd.read_csv(dates_filename, header=0, sep=",")  # raw data
Pop = pd.read_csv(pop_filename, header=0, sep=",")  # raw data
Qr = pd.read_csv(qr_filename, header=0, sep=",")  # raw data
n = Raw.shape[0]
StartFit = Keydates["Dates"][0]
EndFit = Keydates["Dates"][Keydates.shape[0] - 1]
N = Pop["PoblacionMid2020"][np.where(Pop["Entidad"] == entity)[0][0]]
#
di = np.where(Raw["Dates"] == StartFit)[0][0]
df = np.where(Raw["Dates"] == EndFit)[0][0]
A = Raw.values[di:(df + 1), 2].astype(float)  # weekly symptomatic group
H = Raw.values[di:(df + 1), 3].astype(float)  # weekly reported group
D = Raw.values[di:(df + 1), 4].astype(float)  # weekly deaths group
t = np.arange(0, df - di + 1, 1.0)            # observed times
nt = t.shape[0]
#
# Parameter details
#
tcp = Keydates.values[:, 1]
#
tqr = Qr.values[:, 1]
vqr = Qr.values[:, 2]
qr = pchip(tqr, vqr)
#
#
#
ncp = tcp.shape[0] - 1
npar1 = 5  # Number of parameters non-related to change points
npar = npar1 + ncp  # Total number of parameters
#
parnames = ["E0", "I0", "Y0", "alpha", "kappaH"]
for i in range(ncp):
    parnames = parnames + ["b" + str(i + 1)]

Dates0 = pd.date_range(start=StartFit, end=EndFit, freq="1D")
# Fixed parameters
eps = 0.85  # % of symptomatic recovered
rho = 0.2  # % of asymptomatic
q = 0.45  # asymtomtic infection reduction
delta = 1.0 / 7  # asymptomatic recovery rate
eta = 1 / 14  # symptomatic recovery rate
gamma = 1 / 5.1  # incubation rate
kappaA = 1 / 10  # ambulatory recovery rate
# kappaH = 1/10      # hospitalized recovery rate
omega = 1 / 180  # natural immunity
phi = 0  # vaccine imunity
psi = 0  # vaccination rate
chi = 1 / (75.5 * 365)  # Natural death rate
# State specific parameters
State_pars = pd.read_csv(data_path + "State_parameters.csv", header=0, sep=",")
state_index = np.where(State_pars["Entities"] == entity)[0][0]
alpha = State_pars["Alpha"][state_index]  # % of reported recovered
nu = 1.0 / State_pars["Nu"][state_index]  # screening rate
mu = 1.0 / State_pars["Mu"][state_index]  # covid death rate
# qr    = State_pars["qr"][state_index]        # % ambulatory people
if period == "weekly":
    Dates0 = pd.date_range(start=StartFit, end=EndFit, freq="7D")
    delta = delta * 7
    eta = eta * 7
    gamma = gamma * 7
    kappaA = kappaA * 7
    omega = omega * 7
    phi = phi * 7
    psi = psi * 7
    chi = chi * 7
#    nu    = nu*7
#    mu    = mu*7




def ode6(X, t, beta, alpha, kappaH, qr, chi=chi, delta=delta, eps=eps, \
         eta=eta, gamma=gamma, kappaA=kappaA, mu=mu, nu=nu, omega=omega, \
         phi=phi, psi=psi, q=q, rho=rho):
    """ This version calculates cumulative cases for all three classes of
        infected: asymptomatic, symptomatic and reported.
    """
    S = X[0]
    E = X[1]
    I = X[2]
    Y = X[3]
    A = X[4]
    H = X[5]
    R = X[6]
    D = X[7]
    N = S + E + I + Y + R
    N1 = N + A + H
    
    b = beta(t) * (Y + q * I) * S / N
    c1 = qr(t) * nu * (1 - eps) * Y
    c2 = (1 - qr(t)) * nu * (1 - eps) * Y
    c3 = gamma * E * (1 - rho)
    c4 = gamma * E * rho
    dS = omega * R - b + chi * N1 - chi * S
    dE = b - gamma * E - chi * E
    dI = gamma * E * rho - delta * I - chi * I
    dY = gamma * E * (1 - rho) - (eta * eps + nu * (1 - eps)) * Y - \
         chi * Y
    dA = c1 - kappaA * A - chi * A
    dH = c2 - (kappaH * alpha + mu * (1 - alpha)) * H - chi * H
    dR = delta * I + eta * eps * Y + kappaA * A + kappaH * alpha * H - \
         omega * R - chi * R
    dD = mu * (1 - alpha) * H
    CA = c1
    CH = c2
    CY = c3
    CI = c4
    return np.array([dS, dE, dI, dY, dA, dH, dR, dD, CA, CH, CY, CI])





# np.random.seed(4)
MCMC = pytwalk.pytwalk(n=npar, U=logpost, Supp=support)
MCMC.Run(T=ni, x0=rinit(), xp0=rinit())

# MAP
MAP = getMAP(MCMC.Output)

# Energy plot
## energyplot(MCMC.Output[:, -1])

# IAT(MCMC.Output[burnin:ni, :], cols=17)

# Subsampling
Par = getsample(MCMC.Output[:, 0:npar], burnin, ess)
Samples = pd.DataFrame(
    MCMC.Output[
        np.linspace(
            burnin,
            ni,
            num=ess,
            endpoint=False,
            dtype=int), :])
Samples = \
    Samples.append(pd.DataFrame(MAP.reshape(1, npar)))
Samples.to_csv(samples_filename, index=False)

## ess= 100
## Par = getsample(MCMC.Output[:, 0:npar], burnin, ess)
## t raceplot(Par, parnames)
## marginalplot(Par, parnames)
# Probability intervals
## parsummary(Par, parnames)

## R_0 = Par[:, 5]*(rho*gamma/((delta + chi)*(gamma+chi)) + \
##    q*(1-rho)*gamma/((eps*eta + (1-eps)*nu + chi)*(gamma+chi)))

# plt.figure()
# plt.hist(R_0)
# plt.xlabel("Basic reproduction number")

## MAP solution
## betaMAP = interp1d(tcp[1:(ncp+1)], MAP[npar1:npar],\
##
# kind = "linear", bounds_error=False, fill_value = (MAP[npar1], MAP[npar-1]))
# #betaMAP=pchip(tcp, MAP[npar1:npar])
# betaMAP=pchip(tcp, np.concatenate((MAP[npar1], MAP[npar1:npar]), axis=None))
# S0MAP = N - MAP[0]-MAP[1]-MAP[2] 
# X0MAP = np.array([S0MAP, MAP[0], MAP[1], MAP[2], 0, 0, 0, 0, 0, 0, 0, 0])
# MAPsol = odeint(ode6, X0MAP, t, args=(betaMAP, MAP[3]))  ###
# MAPD = np.diff(MAPsol[:, 7])
# MAPA = np.diff(MAPsol[:, 8])
# MAPH = np.diff(MAPsol[:, 9])
# MAPY = np.diff(MAPsol[:, 10])
# MAPI = np.diff(MAPsol[:, 11])

# # Posterior trayectories

# #Dates1 = pd.date_range(start= StartPred, end = EndPred, freq = "7D")
# #Dates2 = pd.date_range(start= StartFit, end = EndPred, freq = "7D")
# #nt1 = Dates2.shape[0]
# #t1 = np.arange(0, nt1, 1)

# # Solutions from posterior samples
# CD = np.zeros(nt*ess).reshape(ess, nt)
# CT = np.zeros(nt*ess).reshape(ess, nt)
# CY = np.zeros(nt*ess).reshape(ess, nt)
# CI = np.zeros(nt*ess).reshape(ess, nt)

# for i in np.arange(0, ess,1):
##beta = interp1d(tcp[1:(ncp+1)], Par[i, npar1:npar],\
##  kind = "linear",
#   bounds_error=False,
#   fill_value = (Par[i, npar1], Par[i, npar-1]))
## beta = pchip(tcp, Par[i, npar1:npar])
# beta =
#   pchip(
#       tcp,
#       np.concatenate((Par[i, npar1], Par[i, npar1:npar]),
#   axis=None))
#   S0 = N - Par[i,0]-Par[i,1]-Par[i,2]
#   X0 = np.array(
#       [S0, Par[i,0], Par[i,1], Par[i,2], 0, 0, 0, 0, 0, 0]
# )
#   sol = odeint(ode6, X0, t, args=(beta, Par[i, 3]))  ###
#   CD[i,:] = sol[:, 6]
#   CT[i,:] = sol[:, 7]
#   CY[i,:] = sol[:, 8]
#   CI[i,:] = sol[:, 9]

# InciD = np.diff(CD[:,], axis=1)
# InciT = np.diff(CT[:,], axis=1)
# InciY = np.diff(CY[:,], axis=1)
# InciI = np.diff(CI[:,], axis=1)


# CD_int = np.quantile(CD, np.array([0.05,0.5,0.95]), axis=0)
# CT_int = np.quantile(CT, np.array([0.05,0.5,0.95]), axis=0)
# CY_int = np.quantile(CY, np.array([0.05,0.5,0.95]), axis=0)
# CI_int = np.quantile(CI, np.array([0.05,0.5,0.95]), axis=0)
# Ctot = CI+CY
# Ctot_int = np.quantile(Ctot, np.array([0.05,0.5,0.95]), axis=0)

# InciD_int  = np.quantile(InciD, np.array([0.05,0.5,0.95]), axis=0)
# InciT_int  = np.quantile(InciT, np.array([0.05,0.5,0.95]), axis=0)
# InciY_int  = np.quantile(InciY, np.array([0.05,0.5,0.95]), axis=0)
# InciI_int  = np.quantile(InciI, np.array([0.05,0.5,0.95]), axis=0)

# plt.figure()
# # Asymptomatic
# plt.subplot(2,2,1)
# plt.plot(Dates0[1:nt], MAPI, "-", linewidth = 3, color="red")
# plt.plot(Dates0[1:nt], InciI_int[0,:], "--", linewidth = 3, color="gray")
# plt.plot(Dates0[1:nt], InciI_int[1,:], "-", linewidth = 3, color="black")
# plt.plot(Dates0[1:nt], InciI_int[2,:], "--", linewidth = 3, color="gray")
# plt.xticks(rotation=-45, ha="left")
# plt.ylabel("Asymptomatic")
# ax=plt.gca()
# ax.xaxis.set_major_locator(ticker.MultipleLocator(15))
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
# plt.grid(True)

# # Symptomatic
# plt.subplot(2,2,2)
# plt.plot(Dates0[1:nt], MAPY, "-", linewidth = 3, color="red")
# plt.plot(Dates0[1:nt], InciY_int[0,:], "--", linewidth = 3, color="gray")
# plt.plot(Dates0[1:nt], InciY_int[1,:], "-", linewidth = 3, color="black")
# plt.plot(Dates0[1:nt], InciY_int[2,:], "--", linewidth = 3, color="gray")
# plt.bar(Dates0, Y, color="tab:blue", width=5)
# plt.xticks(rotation=-45, ha="left")
# plt.ylabel("Symptomatic")
# ax=plt.gca()
# ax.xaxis.set_major_locator(ticker.MultipleLocator(15))
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
# plt.grid(True)

# # Reported
# plt.subplot(2,2,3)
# plt.plot(Dates0[1:nt], MAPT, "-", linewidth = 3, color="red")
# plt.plot(Dates0[1:nt], InciT_int[0,:], "--", linewidth = 3, color="gray")
# plt.plot(Dates0[1:nt], InciT_int[1,:], "-", linewidth = 3, color="black")
# plt.plot(Dates0[1:nt], InciT_int[2,:], "--", linewidth = 3, color="gray")
# plt.bar(Dates0, T, color="tab:blue", width=5)
# plt.xticks(rotation=-45, ha="left")
# plt.ylabel("Reported")
# ax=plt.gca()
# ax.xaxis.set_major_locator(ticker.MultipleLocator(15))
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
# plt.grid(True)

# # Deaths
# plt.subplot(2,2,4)
# plt.plot(Dates0[1:nt], MAPD, "-", linewidth = 3, color="red")
# plt.plot(Dates0[1:nt], InciD_int[0,:], "--", linewidth = 3, color="gray")
# plt.plot(Dates0[1:nt], InciD_int[1,:], "-", linewidth = 3, color="black")
# plt.plot(Dates0[1:nt], InciD_int[2,:], "--", linewidth = 3, color="gray")
# plt.bar(Dates0, D, color="tab:blue", width=5)
# plt.xticks(rotation=-45, ha="left")
# plt.ylabel("Death")
# ax=plt.gca()
# ax.xaxis.set_major_locator(ticker.MultipleLocator(15))
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
# plt.grid(True)

# # Cumulative cases
# plt.figure()
# plt.plot(Dates0, CI_int[0,:], "--", linewidth = 2, color="black")
# plt.plot(Dates0, CI_int[1,:], "-", linewidth = 3, color="black", label = "Asymptomatic")
# plt.plot(Dates0, CI_int[2,:], "--", linewidth = 2, color="black")
# plt.plot(Dates0, CY_int[0,:], "--", linewidth = 2, color="red")
# plt.plot(Dates0, CY_int[1,:], "-", linewidth = 3, color="red", label = "Symptomatic")
# plt.plot(Dates0, CY_int[2,:], "--", linewidth = 2, color="red")
# plt.plot(Dates0, CT_int[0,:], "--", linewidth = 2, color="blue")
# plt.plot(Dates0, CT_int[1,:], "-", linewidth = 3, color="blue", label = "Reported")
# plt.plot(Dates0, CT_int[2,:], "--", linewidth = 2, color="blue")
# plt.plot(Dates0, Ctot_int[0,:], "--", linewidth = 2, color="gray")
# plt.plot(Dates0, Ctot_int[1,:], "-", linewidth = 3, color="gray", label = "Total")
# plt.plot(Dates0, Ctot_int[2,:], "--", linewidth = 2, color="gray")
# plt.xticks(rotation=-45, ha="left")
# plt.ylabel("Cumulative infections")
# ax=plt.gca()
# plt.legend()
# ax.xaxis.set_major_locator(ticker.MultipleLocator(15))
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
# plt.grid(True)


# # Effective contact rates
# if period=="daily":
#     #t = np.arange(0, nt, 0.01)
#     #nt = t.shape[0]
#     Beta = np.zeros((ess, nt))
#     Beta_ch = np.zeros((ess, nt-1))
#     # Contact rate
#     for i in np.arange(0, ess,1):
#         #beta = interp1d(tcp[1:(ncp+1)], Par[i, npar1:npar],\
#         #            kind = "linear", bounds_error=False, fill_value = (Par[i, npar1], Par[i, npar-1]))
#         #beta=pchip(tcp, Par[i, npar1:npar])
#         beta=pchip(tcp, np.concatenate((Par[i, npar1], Par[i, npar1:npar]), axis=None))
#         Beta[i, :] = beta(t)
#         #Beta_ch[i, :] = np.diff(Beta[i,:])

#     Beta_int = np.quantile(Beta, np.array([0.05,0.5,0.95]), axis=0)
#     #Beta_ch_int = np.quantile(Beta_ch, np.array([0.05,0.5,0.95]), axis=0)

#     Beta_int1 = np.quantile(Par[:, npar1:npar], np.array([0.05,0.5,0.95]), axis=0)
#     #Beta_ch_int1 = np.quantile(np.diff(Par[:, npar1:npar]), np.array([0.05,0.5,0.95]), axis=0)

#     plt.figure()
#     #plt.subplot(2,1,1)
#     plt.plot(t, Beta_int[0,:], "--", linewidth = 2, color="black")
#     plt.plot(t, Beta_int[1,:], "-", linewidth = 3, color="black")
#     plt.plot(t, Beta_int[2,:], "--", linewidth = 2, color="black")
#     #plt.plot(tcp[1:(ncp+1)], Beta_int1[0,:], "-", linewidth = 3, color="blue")
#     #plt.plot(tcp[1:(ncp+1)], Beta_int1[1,:], "-", linewidth = 3, color="blue")
#    #plt.plot(tcp[1:(ncp+1)], Beta_int1[2,:], "-", linewidth = 3, color="blue")

#     # plt.plot(Dates0, betaMAP(t), "--", linewidth = 3, color="red")
#     # plt.xticks(rotation=-45, ha="left")
#     # plt.ylabel("Contact rate")
#     # ax=plt.gca()
#     # ax.xaxis.set_major_locator(ticker.MultipleLocator(15))
#     # ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
#     # plt.grid(True)

#     # plt.subplot(2,1,2)
#     # plt.plot(Dates0[1:nt], Beta_ch_int[0,:], "--", linewidth = 2, color="black")
#     # plt.plot(Dates0[1:nt], Beta_ch_int[1,:], "-", linewidth = 3, color="black")
#     # plt.plot(Dates0[1:nt], Beta_ch_int[2,:], "--", linewidth = 2, color="black")
#     # plt.plot(Dates0[1:nt], np.diff(betaMAP(t)), "--", linewidth = 3, color="red")
#     # plt.xticks(rotation=-45, ha="left")
#     # plt.ylabel("Contact rate")
#     # ax=plt.gca()
#     # ax.xaxis.set_major_locator(ticker.MultipleLocator(15))
#     # ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
#     # plt.grid(True)


# #plt.savefig("par.traces",  bbox_inches = 'tight', format="eps")


# # # Optimization solution
# # bounds = [(0,100), (0, 100), (0, 100)]
# # for i in range(1, npar2+1, 1):
# #     bounds = bounds+[(0,1)]
# # MAP1 = differential_evolution(loglike, bounds, workers=-1).x
# # p1 = interp1d(tcp, np.concatenate((MAP1[npar1], MAP1[npar1:(npar1+ncp)], MAP1[npar1+ncp-1]),
# #                                   axis=None))
# # p2 = interp1d(tcp, np.concatenate((MAP1[npar1+ncp], MAP1[(npar1+ncp):(npar1+2*ncp)], MAP1[npar1+2*ncp-1]),
# #                                   axis=None))
# # p3 = interp1d(tcp, np.concatenate((MAP1[npar1+2*ncp], MAP1[(npar1+2*ncp):(npar1+3*ncp)], MAP1[npar1+3*ncp-1]),
# #                                   axis=None))
# # X0MAP1 = np.array([S10, 0, 0, 0, 0, 0, 0, 0, \
# #                    S20, 0, MAP1[0], MAP1[1], MAP1[2], 0, 0, 0, \
# #                    S30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
# # MAPsol1 = odeint(ode, X0MAP1, t, args=(p1,p2,p3))

# # MAP1D1  = np.diff(MAPsol1[:, 7])
# # MAP1D2  = np.diff(MAPsol1[:, 15])
# # MAP1D3  = np.diff(MAPsol1[:, 23])
# # MAP1T1  = np.diff(MAPsol1[:, 24])
# # MAP1T2  = np.diff(MAPsol1[:, 25])
# # MAP1T3  = np.diff(MAPsol1[:, 26])
