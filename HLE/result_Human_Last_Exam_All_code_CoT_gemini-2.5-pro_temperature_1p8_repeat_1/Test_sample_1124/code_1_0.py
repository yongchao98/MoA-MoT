import numpy as np
from scipy.optimize import fsolve

# The task is to modify the logistic map X_n+1 = R * X_n * (1 - X_n)
# to have an equilibrium point near 1.05 for R = 3.57, using only X and R.

# Our proposed modification is: X_n+1 = R * X_n * (1 - X_n) + (R - X_n) / 2

# Set the value for R
R = 3.57

# Define the modified map function f(X)
def modified_logistic_map(x, r_val):
    """
    Calculates the next value (X_n+1) using the modified logistic map.
    """
    return r_val * x * (1 - x) + (r_val - x) / 2

# An equilibrium point X* satisfies the equation f(X*) = X*.
# We can find this by solving for the root of the equation f(X*) - X* = 0.
def find_equilibrium(x, r_val):
    """
    Function to find the root of, representing the equilibrium condition.
    """
    return modified_logistic_map(x, r_val) - x

# We use a numerical solver to find the root of find_equilibrium.
# We provide an initial guess close to our target value of 1.05.
initial_guess = 1.05
# The fsolve function returns an array, so we take the first element.
equilibrium_point = fsolve(find_equilibrium, initial_guess, args=(R))[0]

# --- Output the results ---
print("The modified logistic map is:")
print("X_n+1 = R * X_n * (1 - X_n) + (R - X_n) / 2")
print("-" * 40)
print(f"For R = {R}, the new equilibrium point is found at:")
print(f"X = {equilibrium_point:.6f}")
print("-" * 40)
print("Verification of the equilibrium point:")
print("We check if X_n+1 = X_n for this value.")

# Assign the found equilibrium point to a variable for clarity in the final equation
x_star = equilibrium_point
lhs = x_star
# Calculate the right-hand side of the equation using the equilibrium point
rhs = R * x_star * (1 - x_star) + (R - x_star) / 2

# Print each number in the final equation as requested
print(f"X_n+1 = {R} * {x_star:.6f} * (1 - {x_star:.6f}) + ({R} - {x_star:.6f}) / 2")
print(f"X_n+1 = {rhs:.6f}")
print(f"X_n   = {lhs:.6f}")
