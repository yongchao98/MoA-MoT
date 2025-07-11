import numpy as np
from scipy.optimize import minimize

# Define a strictly convex function f(x) = e^x
# A function is strictly convex if its second derivative is always positive.
# The second derivative of e^x is e^x, which is > 0 for all real x.
def strictly_convex_function(x):
    """A simple strictly convex function."""
    return np.exp(x[0])

# An initial guess for the optimization
x0 = [10.0]

# We will try to find the minimum of this function.
# Because the function f(x) = e^x decreases as x -> -infinity,
# there is no single point 'x' that gives a minimum value.
# The optimization algorithm should fail to converge to a specific minimum.
result = minimize(strictly_convex_function, x0, method='BFGS')

print("--- Attempting to find the minimum of f(x) = e^x ---")
print(f"Was the optimization successful? {result.success}")
print(f"Optimizer message: {result.message}")
print(f"Final 'x' value: {result.x[0]}")
print(f"Final function value (approaches 0): {result.fun}")
print("\nThis demonstrates that a strictly convex function does not necessarily have a global minimizer.")
print("The function f(x) = e^x is a counterexample. The optimization process follows the function's gradient towards negative infinity without ever finding a minimum point.")
print("Therefore, the statement 'Any strictly convex function has a unique global minimizer' is false.")
