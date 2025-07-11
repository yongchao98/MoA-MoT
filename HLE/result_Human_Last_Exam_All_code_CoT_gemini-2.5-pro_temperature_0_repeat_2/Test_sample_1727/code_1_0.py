import numpy as np

def calculate_chi(alpha):
    """
    Calculates the amplitude factor chi as a function of the dimensionless wavenumber alpha = kh.
    
    The full expression is:
    chi = ((2 + alpha) * np.cosh(alpha) - 2 * alpha * np.sinh(alpha)) / (2 * (np.cosh(alpha)**2 + alpha**2))
    
    However, the problem asks for a single stable amplitude factor. This is found by taking the
    long-wavelength limit (alpha -> 0), which is relevant for large-scale tectonic plates.
    """
    
    # In the limit as alpha approaches 0, the expression simplifies to 1.
    # We can demonstrate this by evaluating for a very small alpha.
    if alpha == 0:
        chi = 1.0
    else:
        C = np.cosh(alpha)
        S = np.sinh(alpha)
        numerator = (2 + alpha) * C - 2 * alpha * S
        denominator = 2 * (C**2 + alpha**2)
        chi = numerator / denominator
    return chi

# The problem implies a single value for the factor, which corresponds to the
# long-wavelength limit (alpha -> 0).
chi_stable = calculate_chi(0)

# The final equation is e_s = chi * (S_0 / (Delta_rho * g))
# The problem asks to output the number in the final equation, which is the value of chi.
print(int(chi_stable))