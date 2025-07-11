import numpy as np
from scipy.optimize import fsolve

# Define parameters based on the problem statement.
# The factor beta = 1/(kB*T) is absorbed into the parameters,
# making them dimensionless.
beta_mu = 0.1
beta_epsilon = -1 / (2 * np.pi)
z_horizontal = 4
z_vertical = 8
z = z_horizontal + z_vertical

# The self-consistency equation from mean-field theory for the average
# occupancy 'n' is:
# n = 1 / (1 + exp(beta * (z*epsilon*n - mu)))
# In dimensionless form: n = 1 / (1 + exp(z*beta_epsilon*n - beta_mu))

# We define a function g(n) = 0 that we can solve numerically.
# g(n) = n - (1 / (1 + exp(z*beta_epsilon*n - beta_mu)))
def equation_to_solve(n, z_val, beta_eps, beta_mu_val):
    """Represents the equation g(n) = n - f(n) = 0"""
    return n - 1 / (1 + np.exp(z_val * beta_eps * n - beta_mu_val))

# Provide an initial guess for the average occupancy 'n'.
initial_guess = 0.5

# Use fsolve to find the root of the equation.
# fsolve returns an array, so we take the first element.
solution = fsolve(equation_to_solve, initial_guess, args=(z, beta_epsilon, beta_mu))[0]

# Output the equation with the specific numerical values and the final result.
print("The self-consistency equation to solve is of the form: n = 1 / (1 + exp(z * beta_epsilon * n - beta_mu))")
print("\nUsing the given parameters:")
print(f"  z (total coordination number) = {z}")
print(f"  beta_epsilon (dimensionless interaction energy) = {beta_epsilon:.5f}")
print(f"  beta_mu (dimensionless chemical potential) = {beta_mu}")

print("\nThe final equation with numerical values is:")
print(f"  n = 1 / (1 + exp({z} * ({beta_epsilon:.5f}) * n - {beta_mu}))")

# Print the final calculated answer, rounded to three decimal places.
print(f"\nThe average occupancy per site <n> is: {solution:.3f}")
