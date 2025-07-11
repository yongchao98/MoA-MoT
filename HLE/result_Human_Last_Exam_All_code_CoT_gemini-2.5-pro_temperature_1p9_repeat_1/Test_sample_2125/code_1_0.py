import numpy as np
from scipy import optimize

def equation_for_alpha(alpha):
  """
  Defines the polynomial equation 2*alpha^4 - 2*alpha^3 - 1 = 0.
  The root of this equation gives the value of alpha_0.
  """
  return 2 * alpha**4 - 2 * alpha**3 - 1

# Find the positive root of the equation.
# We can check that f(1) = -1 and f(2) = 15, so a root exists between 1 and 2.
try:
    solution = optimize.root_scalar(equation_for_alpha, bracket=[1, 2])
    alpha_0 = solution.root

    print("The problem reduces to solving the algebraic equation for alpha:")
    print("2 * alpha^4 - 2 * alpha^3 - 1 = 0")
    print("\nSolving this equation numerically gives the value for alpha_0.")
    print(f"The solution is alpha_0 = {alpha_0}")

    # Verification of the solution
    verification_value = 2 * alpha_0**4 - 2 * alpha_0**3 - 1
    print("\nTo verify the solution, we plug it back into the equation:")
    print(f"2 * ({alpha_0:.8f})^4 - 2 * ({alpha_0:.8f})^3 - 1 = {verification_value:.2e}")
except (ImportError, ModuleNotFoundError):
    print("SciPy is not installed. Please install it using 'pip install scipy'.")
    # For a simplified environment without scipy, let's provide a basic implementation
    print("\nAttempting a basic numerical solver (bisection method):")
    a, b = 1.0, 2.0
    tol = 1e-9
    while (b - a) / 2.0 > tol:
        midpoint = (a + b) / 2.0
        if equation_for_alpha(midpoint) == 0:
            a = b = midpoint
            break
        elif np.sign(equation_for_alpha(midpoint)) == np.sign(equation_for_alpha(a)):
            a = midpoint
        else:
            b = midpoint
    alpha_0_manual = (a + b) / 2
    print(f"The solution is alpha_0 = {alpha_0_manual}")

