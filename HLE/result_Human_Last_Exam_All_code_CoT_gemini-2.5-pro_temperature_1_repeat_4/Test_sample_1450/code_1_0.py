import numpy as np
from scipy.optimize import fsolve
from scipy.special import polylog

# --- Explanation ---
print("This script solves for the fugacity 'z' of an ideal Fermi gas.")
print("The condition is that its number density is 75% that of a classical ideal gas at the same pressure P and temperature T.")
print("This leads to a specific equation relating the fugacity 'z' to the ratio of two Fermi-Dirac integrals, f_n(z).\n")

# --- Equation Setup ---
# The number density of a Fermi gas (ν_f) is related to the classical density (ν_c) by:
# ν_f = ν_c * (f_{3/2}(z) / f_{5/2}(z))
# Given ν_f = 0.75 * ν_c, we must solve the equation:
# f_{3/2}(z) / f_{5/2}(z) = 0.75
# Or, f_{3/2}(z) - 0.75 * f_{5/2}(z) = 0

print("The equation to solve for the fugacity 'z' is:")
print("f_3/2(z) / f_5/2(z) = 0.75\n")


# --- Numerical Solution ---

def equation_to_solve(z):
    """
    Defines the equation f_{3/2}(z) - 0.75 * f_{5/2}(z) = 0 for the root finder.
    The Fermi-Dirac integral f_n(z) is related to the polylogarithm function Li_n(x)
    by f_n(z) = -Li_n(-z). We use scipy.special.polylog(n, x) which computes Li_n(x).
    """
    # Ensure z is a float for polylog function
    z_float = float(z)
    
    # Calculate the Fermi-Dirac integrals
    f_3_2 = -polylog(1.5, -z_float)
    f_5_2 = -polylog(2.5, -z_float)
    
    # Return the value of the function whose root we are seeking
    return f_3_2 - 0.75 * f_5_2

# An initial guess for the root. Based on series expansions, the root is around 1.4.
initial_guess = 1.5

# Use fsolve to find the root of the equation
z_solution, = fsolve(equation_to_solve, initial_guess)

# --- Output the Result ---
print(f"The calculated value of the fugacity 'z' that solves the equation is approximately {z_solution:.4f}.")
print("Rounding the result to two significant digits, we get:")
print(f"{z_solution:.1f}")
