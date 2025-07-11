import numpy as np
from scipy.integrate import quad

def estimate_fluctuation_scaling():
    """
    Estimates the scaling function for R(epsilon) by calculating the
    proportionality constant C.
    """

    # The problem is rescaled with X = epsilon*x, so the domain is [0, 1].
    # The leading order fluctuation term Y1(X) is a sum of the scaled
    # Green's function g(X, Z_i). The scaling of R(epsilon) depends on
    # the variance of g(X, Z) where Z is a single random point from U[0,1].

    # Green's function for the scaled operator L[g] = g'' - g' = delta(X-S)
    # on [0,1] with boundary conditions g(0)=0, g(1)=0.
    def greens_g_scaled(x, s):
        e = np.e
        ex = np.exp(x)
        
        if x < s:
            # For x < s, g(x,s) = c2 * (e^x - 1)
            c2 = (1 - np.exp(1 - s)) / (e - 1)
            return c2 * (ex - 1)
        else:
            # For x > s, g(x,s) = d2 * (e^x - e)
            d2 = (1 - np.exp(-s)) / (e - 1)
            return d2 * (ex - e)

    # The expected value E[g(X,Z)] for Z ~ U[0,1] is found by integrating
    # g(X,S) with respect to S from 0 to 1. This can be solved analytically.
    def expected_g_scaled(x):
        e = np.e
        ex = np.exp(x)
        return (ex - 1) / (e - 1) - x

    # The variance is Var[g(X,Z)] = E[g(X,Z)^2] - (E[g(X,Z)])^2.
    # We compute E[g(X,Z)^2] by numerically integrating g(X,S)^2 over S.
    def variance_g_scaled(x):
        # E[g(X,Z)^2] = integral from 0 to 1 of g(X,S)^2 dS
        integrand_sq = lambda s: greens_g_scaled(x, s)**2
        E_g_sq, _ = quad(integrand_sq, 0, 1, limit=200)
        
        # E[g(X,Z)] from the analytical formula
        E_g = expected_g_scaled(x)
        
        return E_g_sq - E_g**2

    # To find the constant C, we need C^2 = max_X Var[g(X,Z)].
    # We search for this maximum over a grid of X values in [0, 1].
    x_grid = np.linspace(0.01, 0.99, 200)
    var_values = np.array([variance_g_scaled(x) for x in x_grid])

    # The maximum variance gives C^2.
    max_var = np.max(var_values)
    C = np.sqrt(max_var)
    
    # The final result is the function R(epsilon)
    scaling_exponent = 0.5
    print("The estimated relationship for the maximum magnitude of fluctuations is:")
    print(f"R(epsilon) = C * epsilon^k")
    print(f"where the scaling exponent k = {scaling_exponent}")
    print(f"and the constant C is sqrt(max_X Var[g(X,Z)])")
    print(f"Numerically, C ~= {C:.4f}")
    print("\nFinal equation:")
    print(f"R(epsilon) ~= {C:.4f} * epsilon^{scaling_exponent}")

estimate_fluctuation_scaling()