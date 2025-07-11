import numpy as np
from scipy.optimize import root_scalar
from scipy.special import polylog

def solve_fugacity():
    """
    This function solves for the fugacity z of a Fermi gas under the specified conditions.
    """
    # The problem reduces to solving the equation f_3/2(z) = 0.75 * f_5/2(z).
    # We define a function whose root we want to find.
    def equation_to_solve(z):
        """
        Represents the equation f_3/2(z) - 0.75 * f_5/2(z) = 0.
        The Fermi-Dirac integral f_n(z) is computed as -polylog(n, -z).
        """
        # Fermi-Dirac integral f_n(z)
        f_3_2 = -polylog(1.5, -z)
        f_5_2 = -polylog(2.5, -z)
        
        # We only care about the real part of the solution.
        return np.real(f_3_2 - 0.75 * f_5_2)

    try:
        # Use a numerical solver to find the root of the equation.
        # An approximation suggests the root is near 1.1, so we search in a nearby bracket.
        solution = root_scalar(equation_to_solve, bracket=[0.5, 1.5], method='brentq')
        z_value = solution.root

        # Print the derivation and the final equation with its numbers.
        print("The number density of the Fermi gas (nu_F) is 75% that of a classical gas (nu_cl).")
        print("This condition leads to the following relationship between the Fermi-Dirac integrals (f_n) and the fugacity (z):")
        print("f_3/2(z) = 0.75 * f_5/2(z)\n")
        
        print("Solving this equation numerically gives the value of the fugacity z.")
        print("\nThe value of the fugacity is:")
        # Format the result to two significant digits using the '.2g' format specifier.
        print(f"{z_value:.2g}")

    except ImportError:
        print("Error: The 'scipy' library is required. Please install it using 'pip install scipy'.")
    except Exception as e:
        print(f"An error occurred during calculation: {e}")

if __name__ == '__main__':
    solve_fugacity()