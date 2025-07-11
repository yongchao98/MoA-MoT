import numpy as np
from scipy.integrate import quad
from scipy.special import gamma
from scipy.optimize import root_scalar

def solve_fermi_fugacity():
    """
    This script calculates the fugacity 'z' for an ideal Fermi gas under the condition
    that its number density is 75% that of a classical ideal gas at the same
    pressure and temperature. This condition simplifies to the equation:
    f_{1.5}(z) = 0.75 * f_{2.5}(z)
    where f_n(z) is the Fermi-Dirac integral of order n.
    """

    # --- Step 1: Define a robust function for the Fermi-Dirac integral ---
    # The integral is calculated numerically for accuracy, especially for z > 1.
    def fermi_dirac_integral(n, z):
        """
        Computes the Fermi-Dirac integral f_n(z) using numerical integration.
        f_n(z) = (1/Gamma(n)) * integral from 0 to inf of [x^(n-1) / (exp(x)/z + 1)] dx
        """
        if z < 0:
            raise ValueError("Fugacity z must be non-negative.")
        if z == 0:
            return 0.0

        # Define the integrand for the numerical integration.
        # np.exp(x - np.log(z)) is a numerically stable way to write exp(x)/z.
        def integrand(x):
            return (x**(n - 1)) / (np.exp(x - np.log(z)) + 1)

        # Perform the numerical integration from 0 to infinity
        integral_val, _ = quad(integrand, 0, np.inf)

        return integral_val / gamma(n)

    # --- Step 2: Define the equation to be solved for the root 'z' ---
    # We need to find the root of the function: F(z) = f_{1.5}(z) - 0.75 * f_{2.5}(z) = 0
    
    # Store the numbers from the final equation as variables
    n1 = 1.5
    n2 = 2.5
    factor = 0.75

    def equation_to_solve(z):
        """Represents the equation F(z) = 0 that we need to solve."""
        val1 = fermi_dirac_integral(n1, z)
        val2 = fermi_dirac_integral(n2, z)
        return val1 - factor * val2

    # --- Step 3: Solve the equation numerically ---
    # We use a root-finding algorithm to find the value of z.
    # Preliminary analysis shows the root is between 1 and 4.

    # Print the equation we are solving, including each number as requested.
    print(f"Solving the equation for the fugacity 'z': f_{n1}(z) = {factor} * f_{n2}(z)")

    try:
        # Use the root_scalar function to find the root within the bracket [1, 4].
        # Brent's method (brentq) is robust and efficient for a bracketed root.
        solution = root_scalar(equation_to_solve, bracket=[1, 4], method='brentq')

        if solution.converged:
            fugacity_z = solution.root
            # Print the final result formatted to two significant digits.
            print(f"\nThe value of the fugacity 'z' is: {fugacity_z:.2g}")
        else:
            print("\nThe solver did not converge to a solution.")
            print(f"Solver details: {solution}")

    except ValueError as e:
        print(f"An error occurred during calculation: {e}")

if __name__ == '__main__':
    solve_fermi_fugacity()