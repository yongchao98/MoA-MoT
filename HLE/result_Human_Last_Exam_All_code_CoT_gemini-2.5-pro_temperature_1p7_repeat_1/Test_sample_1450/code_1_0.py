import numpy as np
from scipy.integrate import quad
from scipy.special import gamma
from scipy.optimize import root_scalar

def solve_fugacity():
    """
    This script solves for the fugacity z of a Fermi gas under the specified conditions.
    """
    # Step 1: Define the Fermi-Dirac integral function f_n(z) using numerical integration.
    # The integrand is x^(n-1) / (exp(x)/z + 1).
    def integrand(x, n, z):
        # Change of variable to handle the step-like behavior for large z
        # Let alpha = ln(z). Denominator is exp(x-alpha) + 1.
        # This helps the numerical integrator converge.
        if z <= 0:
            return 0.0
        alpha = np.log(z)
        if x > alpha:
             # For x > alpha, use identity to avoid overflow in exp(x-alpha)
            return (x**(n-1)) * np.exp(alpha - x) / (1 + np.exp(alpha - x))
        else:
            return (x**(n-1)) / (np.exp(x - alpha) + 1)

    def fd_integral(n, z):
        """Calculates the Fermi-Dirac integral f_n(z)."""
        if z <= 0:
            return 0.0
        # The integral is from 0 to infinity. quad returns a tuple (result, error).
        result, _ = quad(integrand, 0, np.inf, args=(n, z))
        return result / gamma(n)

    # Step 2: Define the target function whose root (z) we need to find.
    # The condition is nu_fermi = 0.75 * nu_classical
    # which simplifies to f_3/2(z) = 0.75 * f_5/2(z).
    ratio = 0.75
    def root_function(z, n1, n2, target_ratio):
        """The function that should be zero at the solution."""
        return fd_integral(n1, z) - target_ratio * fd_integral(n2, z)

    try:
        # Step 3: Use a root-finding algorithm to solve for z.
        # Based on textbook tables, the root is expected to be around z=3.
        # We search in the bracket [1, 5] where the function changes sign.
        solution = root_scalar(root_function, args=(1.5, 2.5, ratio), bracket=[1, 5])
        z_solution = solution.root
        
        # Step 4: Print the results as requested.
        print("The problem requires solving the following equation for the fugacity z:")
        print(f"f_3/2(z) = {ratio} * f_5/2(z)")
        print("where f_n(z) is the Fermi-Dirac integral of order n.\n")
        
        # Calculate the values for the final equation display
        f32_val = fd_integral(1.5, z_solution)
        f52_val = fd_integral(2.5, z_solution)
        
        # Format the result to two significant digits using the 'g' format specifier
        print(f"The numerical solution for the fugacity is z = {z_solution:.2g}\n")
        
        print("Substituting this value back verifies the equation:")
        print(f"f_3/2({z_solution:.2g}) ≈ {f32_val:.4f}")
        print(f"{ratio} * f_5/2({z_solution:.2g}) ≈ {ratio} * {f52_val:.4f} ≈ {ratio * f52_val:.4f}\n")
        
        # Final output of the equation with numbers
        print("Final Equation:")
        print(f"{f32_val:.4f} = {ratio * f52_val:.4f}")

    except Exception as e:
        print(f"An error occurred during calculation: {e}")

if __name__ == '__main__':
    solve_fugacity()