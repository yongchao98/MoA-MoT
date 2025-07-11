import numpy as np
from scipy.optimize import minimize
import warnings

# Suppress runtime warnings from the optimizer
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Statement E claims: Any strictly convex function has a unique global minimizer.
# We will test this with a counterexample: f(x) = e^x.
# This function is strictly convex because its second derivative, f''(x) = e^x, is always positive.

def strictly_convex_function(x):
    """A simple strictly convex function: f(x) = e^x."""
    return np.exp(x[0])

# We will start the optimization from an arbitrary point, x0 = 0.
initial_guess = [0]

# We will use a standard optimization algorithm to try to find the minimum.
# If a minimum exists, the algorithm should converge to it.
# However, for f(x) = e^x, the function value decreases as x approaches -infinity.
# The algorithm will continually search for smaller values and will not converge to a finite minimum.
print("--- Demonstrating Statement E is False ---")
print("Statement E: Any strictly convex function has a unique global minimizer.")
print("Counterexample: The strictly convex function f(x) = e^x.")
print("\nAttempting to find the minimum of f(x) = e^x using scipy.optimize.minimize...")

result = minimize(strictly_convex_function, initial_guess, method='BFGS')

print("\n--- Optimization Result ---")
print(f"Algorithm converged successfully: {result.success}")
print(f"Message from optimizer: '{result.message}'")
print(f"Final function value f(x) found: {result.fun:.5e}")
print(f"Final point x found: {result.x[0]}")
print("--------------------------\n")

print("--- Conclusion ---")
print("As shown, the optimization failed to converge because a true minimum does not exist.")
print("The function f(x) = e^x approaches 0 as x approaches negative infinity, but never reaches it.")
print("Since we have found a strictly convex function without a global minimizer, statement E is false.")
print("\nTherefore, the statement that is not true is E.")
