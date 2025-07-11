import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

def estimate_fluctuation_scaling():
    """
    This function solves for the scaling of the fluctuation R(epsilon).
    The logic is explained in the text above. Here we compute the constant
    factor in the scaling law R = C * epsilon**(1/2).
    """

    # Euler's number
    e = np.e

    def G(X, s):
        """
        Calculates the Green's function G(X, s) for the operator
        L = d^2/dX^2 - d/dX with zero boundary conditions on [0, 1].
        """
        # The derivation splits the function based on whether s is less than or greater than X.
        if X < s:
            # Formula for G(X, s) when s > X
            return ((1 - e**(1-s)) / (e - 1)) * (e**X - 1)
        else: # X >= s
            # Formula for G(X, s) when s < X
            return ((1 - e**(-s)) / (e - 1)) * (e**X - e)

    def integrand(s, X):
        """The function to be integrated: G(X,s)^2."""
        return G(X, s)**2

    def variance_integral(X):
        """
        Calculates the integral of G(X,s)^2 over s from 0 to 1.
        The integral must be split at the point s = X because the
        definition of G(X,s) changes there.
        """
        if X <= 1e-9 or X >= 1-1e-9:
            return 0
        
        # Integrate from 0 to X
        integral1, _ = quad(integrand, 0, X, args=(X,))
        # Integrate from X to 1
        integral2, _ = quad(integrand, X, 1, args=(X,))
        return integral1 + integral2

    # To find the maximum of variance_integral(X), we can find the minimum
    # of its negative. We search for the maximum in the interval (0, 1).
    res = minimize_scalar(lambda x: -variance_integral(x), bounds=(0, 1), method='bounded')

    # The maximum value of the integral part of the variance
    max_integral_val = -res.fun
    
    # The value C in R = C * epsilon**(1/2)
    C = np.sqrt(max_integral_val)

    # Output the final estimated relationship for R(epsilon)
    print("The maximum magnitude of fluctuations R is estimated to scale with epsilon as:")
    print(f"R(epsilon) = C * epsilon^(1/2)")
    print("where the constant C is calculated below.")
    print(f"C = (max_X integral[0 to 1] G(X,s)^2 ds)^(1/2) ≈ {C:.4f}")
    
    print("\nFor the second case where z_i ~ Normal(i, 0.5):")
    print("Yes, the scaling for R(epsilon) is expected to remain the same, i.e., R ~ epsilon^(1/2).")
    print("This is because the scaling is primarily determined by the number of sources (N ~ 1/epsilon), the strength of each source (epsilon^2), and the properties of the Green's function, which are robust to the change in the source distribution.")

estimate_fluctuation_scaling()
<<<R(epsilon) = 0.1783 * epsilon^(1/2), Yes>>>