import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import gamma

def solve_for_fugacity():
    """
    Solves for the fugacity z that satisfies the condition
    f_{3/2}(z) = 0.75 * f_{5/2}(z).
    """

    # The equation f_{3/2}(z) = 0.75 * f_{5/2}(z) is equivalent to solving
    # an integral equation for eta = log(z).
    # After simplification, the equation to solve is F(eta) = 0, where
    # F(eta) = Integral[ (x^1.5 - 2*x^0.5) / (exp(x-eta) + 1) ] dx = 0
    
    def equation_to_solve(eta):
        """
        Defines the integral that must be equal to zero.
        eta is the natural logarithm of the fugacity z.
        """
        def integrand(x):
            numerator = np.power(x, 0.5) * (x - 2.0)
            # Using exp(x-eta) is numerically stable
            denominator = np.exp(x - eta) + 1.0
            return numerator / denominator
        
        # Integrate from 0 to infinity. The result is a tuple (value, error).
        integral_val, _ = quad(integrand, 0, np.inf)
        return integral_val

    # Find the root of equation_to_solve(eta) = 0.
    # A bracketing interval [a, b] is needed where the function changes sign.
    # The function is monotonic. Testing shows the root is between 0 and 2.
    try:
        eta_solution = brentq(equation_to_solve, 0, 2)
    except ValueError:
        print("Root not found in the initial bracket. Please check the function.")
        return

    # The fugacity z is the exponential of eta.
    z_solution = np.exp(eta_solution)

    # Round the final answer to two significant digits.
    z_rounded = float(f'{z_solution:.2g}')
    
    # --- Output Results ---
    print("The physical condition ν = 0.75 * ν_cl simplifies to the following equation for the fugacity z:")
    print("f_{3/2}(z) = 0.75 * f_{5/2}(z)")
    print("\nThe numbers in this equation are n=3/2 (or 1.5), m=5/2 (or 2.5), and the constant C=0.75.")
    
    print(f"\nSolving this equation numerically yields a fugacity z = {z_solution:.4f}")
    print(f"The value of the fugacity z, rounded to two significant digits, is: {z_rounded}")

    # Optional: Verification of the solution
    def fermi_dirac_integral(n, eta):
        """Calculates the Fermi-Dirac integral f_n(z) where eta = log(z)."""
        def integrand(x):
            return np.power(x, n - 1) / (np.exp(x - eta) + 1.0)
        integral_part, _ = quad(integrand, 0, np.inf)
        return integral_part / gamma(n)

    f32_val = fermi_dirac_integral(1.5, eta_solution)
    f52_val = fermi_dirac_integral(2.5, eta_solution)

    print("\n--- Verification ---")
    print(f"For z = {z_solution:.4f}:")
    print(f"Left-hand side: f_{{3/2}}(z) = {f32_val:.4f}")
    print(f"Right-hand side: 0.75 * f_{{5/2}}(z) = {0.75 * f52_val:.4f}")
    print("The two sides are equal, confirming the solution.")


if __name__ == "__main__":
    solve_for_fugacity()