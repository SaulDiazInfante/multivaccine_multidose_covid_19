npar = 5
npar1 = 6


def support(theta, n=npar, n1=npar1):
    """ Returns whether the vector theta is in the support or not.
        Parameters:
        ----------
            theta[0]   = E10
            theta[1]   = I10
            theta[2]   = Y10
            theta[3]   = alpha
            theta[4]   = kappaH
            theta[5:n] = bi
        Returns:
        -------
            True or False
    """
    return all(theta > 0) and all(theta[n1:n] < 5) and all(theta[0:3] < 10)\
           and theta[3] < 1 and theta[4] < 20


def logprior(theta, a=20, b=2 / 0.9):
    """ Evaluates the joint prior distribution.
        Parameters:
        --------
            theta[0]   = E10
            theta[1]   = I10
            theta[2]   = Y10
            theta[3]   = alpha
            theta[4]   = kappaH
            theta[5:n] = bi
    """
    return (a - 1) * np.log(theta) - (b - 1) * np.log(1 - theta)


def logpost(theta, N=N, qr=qr):
    """ Evaluates the log-posterior function.
        Parameters:
        --------
        theta[0]   = E10
        theta[1]   = I10
        theta[2]   = Y10
        theta[3:n] = bi
    """
    # betas = np.concatenate((theta[npar1],theta[npar1:npar]), axis=None)
    # beta = interp1d(tcp[1:(ncp+1)], theta[npar1:npar],\
    # kind = "linear", bounds_error = False,
    # fill_value = (theta[npar1], theta[npar-1]))
    beta = pchip(
            tcp,
            np.concatenate((theta[npar1], theta[npar1:npar]),axis=None)
        )
    # beta = interp1d(tcp[1:(ncp+1)], theta[npar1:npar],\
    # kind = "linear", bounds_error=False
    # fill_value = (theta[npar1], theta[npar-1]))
    S0 = N - theta[0] - theta[1] - theta[2]
    X0 = np.array([S0, theta[0], theta[1], theta[2], 0, 0, 0, 0, 0, 0])
    sol = odeint(ode4, X0, t, args=(beta, theta[3], 1 / theta[4], qr))
    mu1 = np.diff(sol[:, 7])
    mu2 = np.diff(sol[:, 8])
    mu3 = np.diff(sol[:, 9])
    # return -(D[1:nt]*np.log(mu1) - mu1).sum() - logprior(theta[3])
    res = -(
            D[1:nt] * np.log(mu1)
            + A[1:nt] * np.log(mu2)
            + H[1:nt] * np.log(mu3)
            - mu1 - mu2 - mu3).sum() \
          - logprior(theta[3])
    return res


def rinit(npar1) :
    """ Returns initial values for each parameter.
        Parameters:
        ----------
        theta[0]  = E10
        theta[1]  = I10
        theta[2]  = Y10
        theta[3:n] = bi
    """
    flag = False
    k = 0
    condition = ((flag is False) and (k < 1000))
    while condition is True:
        if k == 1000:
            print("Initial points couldn't be found")
            exit()
        theta = np.zeros(eta)
        theta[0:3] = np.random.uniform(0, 10, size=3)
        theta[3] = np.random.uniform(0, 1, size=1)
        theta[4] = np.random.uniform(1, 20, size=1)
        theta[npar1:eta] = np.random.uniform(0, 5, size=eta - npar1)
        if np.isnan(logpost(theta)):
            k = k + 1
        else:
            flag = True
        condition = ((flag is False) and (k < 1000))
    return theta
