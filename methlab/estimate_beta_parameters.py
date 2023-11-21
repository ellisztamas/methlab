def estimate_beta_parameters(mu, sigma):
    """Estimate shape parameters of a beta distribution.

    Estimate shape parameters a and b of a beta distribution using the mean and
    variance of that distribution by method-of-moments.

    Parameters
    ----------
    mu : float
        Mean of the beta distribution.
    sigma : float
        Variance of the beta distribution.

    Returns
    -------
    tuple
        Tuple with two elements giving shape parameters a and b in that order.
    """
    k = ((mu * (1-mu)) / sigma) -1
    a = mu * k
    b = (1-mu) * k
    return (a,b)