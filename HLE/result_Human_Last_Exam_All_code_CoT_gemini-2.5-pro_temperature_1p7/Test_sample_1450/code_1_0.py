import numpy as np
from scipy.special import polylog
from scipy.optimize import fsolve

def fermi_dirac_integral(n, z):
    """
    Computes the Fermi-Dirac integral f_n(z) using the polylogarithm function.
    f_n(z) = -Li_n(-z)
    """
    # The polylog function can be complex, so we take the real part.
    # For positive z, -z is negative, and Li_n is real.
    return -np.real(polylog(n, -z))

def equation_to_solve(z):
    """
    The equation we need to find the root of:
    f_{3/2}(z) - 0.75 * f_{5/2}(z) = 0
    """
    f_3_2 = fermi_dirac_integral(1.5, z)
    f_5_2 = fermi_dirac_integral(2.5, z)
    return f_3_2 - 0.75 * f_5_2

# An initial guess for the fugacity z. z > 0.
# A guess around 1.0 is reasonable for this non-degenerate to mildly degenerate regime.
initial_guess = 1.0

# Use fsolve to find the root of the equation
# fsolve returns an array, so we take the first element
fugacity_solution = fsolve(equation_to_solve, initial_guess)[0]

# Print the equation that was solved, including all numbers
n1 = 3/2
n2 = 5/2
factor = 0.75
print(f"The equation being solved is: f_{n1}(z) = {factor} * f_{n2}(z)")

# Print the final result formatted to two significant digits
print(f"The calculated value of the fugacity z is: {fugacity_solution:.2g}")

<<<1.1>>>