#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 17:28:27 2021

@author: mario
"""

import numpy as np
import matplotlib.pyplot as plt


def energyplot(x):
    """ Creates a plot of the energy function.
    """
    plt.figure()
    plt.plot(x)
    plt.xlabel(r"$t$")
    plt.ylabel(r"$- \ln(\pi(x_t|y))$")
    plt.title("log-posterior plot")
    plt.tight_layout()
    plt.show(block=False)
    
    
def traceplot(x, parnames):
    """ Marginal traces
    """
    d = x.shape[1]
    flag = False
    nrow = 2
    while flag == False:
        ncol = np.ceil(d/nrow)
        if ncol >6:
            nrow = nrow+1
        else:
            flag=True
    plt.figure()
    for i in range(0, d, 1):
        plt.subplot(nrow,ncol,i+1)
        plt.plot(x[:,i])
        plt.ylabel(parnames[i])
        plt.xlabel("Iteration")
        plt.axhline(y = np.median(x[:,i]), color="red")
        plt.grid()
    

def marginalplot(x, parnames):
    """ Marginal empirical posterior distributions
    """
    d = x.shape[1]
    flag = False
    nrow = 2
    while flag == False:
        ncol = np.ceil(d/nrow)
        if ncol >6:
            nrow = nrow+1
        else:
            flag=True
    plt.figure()
    for i in range(0, d, 1):
        plt.subplot(nrow,ncol,i+1)
        plt.hist(x[:,i], density = True)
        plt.ylabel(parnames[i])
        plt.xlabel("Iteration")
        plt.axvline(x = np.median(x[:,i]), color="red")
        plt.grid()


def getMAP(x):
    """
    """
    lpost = x[:, -1]
    npar = x.shape[1]-1
    MAPindex = np.where(lpost ==lpost.min())
    MAP = x[MAPindex[0][0], 0:npar]
    return MAP


def getsample(x, burnin, ess):
    """
    """
    n = x.shape[0]
    Par = x[np.linspace(burnin, n, num=ess, endpoint=False, dtype=int), :]
    return Par

def parsummary(x, parnames, alpha=0.05):
    """
    """
    d = x.shape[1]
    Par_int = np.quantile(x, np.array([alpha/2,0.5,1-alpha/2]), axis=0)
    for i in range(0, d, 1):
        print(parnames[i]+": ", (Par_int[:, i]).round(3))
    return

def R_0(beta, rho, gamma, delta, chi, q, eps, eta, nu):
    """
    """
    r0 = beta*(rho*gamma/((delta + chi)*(gamma+chi)) + \
    q*(1-rho)*gamma/((eps*eta + (1-eps)*nu + chi)*(gamma+chi)))
    plt.figure()
    plt.hist(r0)
    plt.xlabel("Basic reproduction number")
    return