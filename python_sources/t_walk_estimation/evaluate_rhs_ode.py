"""
    This module computes the Hermite polynomial to
    evaluate the infection rate beta(t) and symptomatic
    proportion q(t)
"""
import numpy as np


def beta(t):
    return 0.001 * t ** 2


def q_r(t):
    return 0.001 * t ** 2


def evaluate_rhs_ode(x, t, **kwargs):
    """ This version calculates cumulative reported cases and it
        does not considers
        vaccination.
        beta = kwargs['beta']
    """
    alpha = kwargs['alpha']
    kappa_h = kwargs['kappa_h']
    chi = kwargs['chi']
    delta = kwargs['delta']
    eps = kwargs['eps']
    eta = kwargs['eta']
    gamma = kwargs['gamma']
    omega = kwargs['omega']
    kappa_a = kwargs['kappa_a']
    mu = kwargs['mu']
    nu = kwargs['nu']
    omega = kwargs['omega']
    phi = kwargs['phi']
    psi = kwargs['psi']
    q = kwargs['q']
    rho = kwargs['rho']  
            
    s = x[0]
    e = x[1]
    i = x[2]
    y = x[3]
    a = x[4]
    h = x[5]
    r = x[6]
    d = x[7]
    n = s + e + i + y + r
    n1 = n + a + h
    
    b = beta(t) * (y + q * i) * s / n
    ds = omega * r - b + chi * n1 - chi * s
    de = b - gamma * e - chi * e
    di = gamma * e * rho - delta * i - chi * i
    dy = gamma * e * (1 - rho) - (eta * eps + nu * (1 - eps)) * y - chi * y
    da = q_r(t) * nu * (1 - eps) * y - kappa_a * a - chi * a
    dh = (1 - q_r(t)) * nu * (1 - eps) * y - \
         (kappa_h * alpha + mu * (1 - alpha)) * h - chi * h
    dr = delta * i + eta * eps * y + kappa_a * a + kappa_h * alpha * h - \
         omega * r - chi * r
    dd = mu * (1 - alpha) * h
    Ca = q_r(t) * nu * (1 - eps) * y
    Ch = (1 - q_r(t)) * nu * (1 - eps) * y
    return np.array([ds, de, di, dy, da, dh, dr, dd, Ca, Ch])
