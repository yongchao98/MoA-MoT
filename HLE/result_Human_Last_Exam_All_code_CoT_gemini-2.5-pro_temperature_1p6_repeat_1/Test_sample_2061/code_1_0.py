import numpy as np
from scipy.optimize import fsolve

# Given parameters from the problem
T_val = np.log(10)
B_val = 0.5 * 10**20 / (99**2)
alpha_param_A = 2.0
alpha_param_K = 3.0

# Calculate intermediate constants based on T
c2 = 1 - np.exp(-2 * T_val)
c3 = 1 - np.exp(-3 * T_val)

# Define the coefficients of the polynomial in alpha.
# Based on the derived relation B = A^8/8 - K*A^4/4 which allows a positive B.
C8 = (alpha_param_A**8 / 8) / c2**8
C5 = (alpha_param_K * alpha_param_A**4 / 4) / (c3 * c2**4)

def equation_to_solve(alpha):
    """
    This function represents the polynomial equation for alpha:
    C8 * alpha^8 - C5 * alpha^5 - B = 0
    """
    return C8 * alpha**8 - C5 * alpha**5 - B_val

# For a numerical solver, a good initial guess is helpful.
# For large alpha, the alpha^8 term dominates, so alpha^8 is approx. B/C8.
initial_guess = (B_val / C8)**(1/8)

# Use a numerical solver to find the root
alpha_solution = fsolve(equation_to_solve, initial_guess)[0]

# Print the final equation with numerical coefficients for clarity
print("The problem reduces to solving the following equation for alpha:")
# Recompute A and K coefficients using given numbers
A_coeff = 2.0
K_coeff = 3.0
term8_coeff = (A_coeff**8 / 8) / ((1 - 10**-2)**8)
term5_coeff = (K_coeff * A_coeff**4 / 4) / ((1 - 10**-3) * (1 - 10**-2)**4)
print(f"{term8_coeff:.4e} * alpha^8 - {term5_coeff:.4e} * alpha^5 = {B_val:.4e}")

# Print the final answer
print("\nThe solution for alpha is:")
print(alpha_solution)