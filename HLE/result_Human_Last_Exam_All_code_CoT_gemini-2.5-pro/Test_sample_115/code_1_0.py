import numpy as np
from scipy.optimize import minimize

# We want to find the growth rate alpha.
# Based on the analysis, alpha is the square root of the maximum eigenvalue of a 2x2 matrix.
# This eigenvalue, lambda, depends on a parameter y = sin^2(theta), where y is in [0, 1].
# The function for the largest eigenvalue is lambda(y) = (1+y + sqrt(1+2y-3y^2))/2.
# We need to find the maximum value of this function for y in [0, 1].

def neg_lambda(y):
    """
    We want to maximize lambda(y), which is equivalent to minimizing -lambda(y).
    This function returns the negative of the largest eigenvalue.
    y corresponds to sin^2(theta).
    """
    if not (0 <= y <= 1):
        return np.inf  # Return a large value for out-of-bounds input
    discriminant = 1 + 2 * y - 3 * y**2
    if discriminant < 0:
        # This should not happen for y in [0,1]
        return np.inf
    
    eigenvalue = (1 + y + np.sqrt(discriminant)) / 2
    return -eigenvalue

# We use a numerical optimizer to find the minimum of -lambda(y), which corresponds to the maximum of lambda(y).
# We provide an initial guess (x0) and bounds for y.
initial_guess = 0.5
bounds = [(0, 1)]

result = minimize(neg_lambda, initial_guess, bounds=bounds)

# The maximum eigenvalue is the negative of the minimum value found by the optimizer.
max_eigenvalue = -result.fun
y_max = result.x[0]

# The growth rate alpha is the square root of the maximum eigenvalue.
alpha = np.sqrt(max_eigenvalue)

print(f"The analysis of the recursive structure leads to a growth factor alpha which is the square root of the maximum of a function lambda(y).")
print(f"The function to maximize is lambda(y) = (1+y + sqrt(1+2y-3y^2))/2 for y in [0, 1].")
print(f"Numerical optimization finds the maximum value of lambda(y) at y = {y_max:.4f}.")
print(f"The maximum value of the eigenvalue is {max_eigenvalue:.4f}.")
print(f"The exact maximum value is 4/3.")
print(f"The growth rate alpha is the square root of this maximum value.")
print(f"alpha = sqrt(4/3) = 2/sqrt(3)")
print(f"The calculated value of alpha is: {alpha}")
