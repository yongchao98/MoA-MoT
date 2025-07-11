import numpy as np
from scipy.optimize import fsolve

# --- Parameters ---
# We are given parameters in terms of k_B*T, so we can work with dimensionless quantities.
# beta = 1 / (k_B * T)
# beta * mu
beta_mu = 0.1
# beta * epsilon
beta_epsilon = -1 / (2 * np.pi)
# Coordination numbers
z_horizontial = 4
z_vertical = 8
# Total coordination number for a bulk site
z = z_horizontial + z_vertical

# --- Self-Consistency Equation ---
# The equation is: <n> = 1 / (1 + exp(-beta * (mu + z * epsilon * <n>)))
# Let n be <n>. The equation to solve is f(n) = 0 where:
# f(n) = n - 1 / (1 + exp(-(beta_mu + z * beta_epsilon * n)))

def equation_to_solve(n):
    """
    Defines the self-consistency equation in the form f(n) = 0.
    n represents the average occupancy <n>.
    """
    exponent = -(beta_mu + z * beta_epsilon * n)
    return n - 1 / (1 + np.exp(exponent))

# --- Solve Numerically ---
# Initial guess for the average occupancy <n>. Must be between 0 and 1.
initial_guess = 0.5
# Use scipy's fsolve to find the root of the equation.
solution = fsolve(equation_to_solve, initial_guess)
# The result is an array, so we take the first element.
avg_occupancy = solution[0]

# --- Output the results ---
# The prompt asks to output the numbers in the final equation.
# Equation is <n> = 1 / (1 + exp(-(c1 + c2 * <n>)))
c1 = beta_mu
c2 = z * beta_epsilon
print("Solving the self-consistency equation for the average occupancy <n>:")
print(f"<n> = 1 / (1 + exp(-({c1:.3f} + ({c2:.3f}) * <n>)))")
print("\nNumerically solving this equation gives:")
print(f"Average occupancy per site <n> = {avg_occupancy:.3f}")
