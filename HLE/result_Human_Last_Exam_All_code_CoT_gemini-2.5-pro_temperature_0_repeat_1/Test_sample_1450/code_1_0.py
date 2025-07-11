import numpy as np
from scipy.special import polylog
from scipy.optimize import fsolve

def find_fugacity():
    """
    Calculates the fugacity z for an ideal Fermi gas under the condition that its
    number density is 75% that of a classical ideal gas at the same P and T.

    The condition simplifies to the equation:
    f_{3/2}(z) = 0.75 * f_{5/2}(z)

    where f_n(z) is the Fermi-Dirac integral. This is solved numerically.
    """

    # The equation to solve for z.
    # We use the relation f_n(z) = -polylog(n, -z).
    # The equation becomes:
    # -polylog(1.5, -z) = 0.75 * (-polylog(2.5, -z))
    # polylog(1.5, -z) - 0.75 * polylog(2.5, -z) = 0
    def equation_to_solve(z):
        if z <= 0:
            return np.inf # Fugacity must be positive
        return polylog(1.5, -z) - 0.75 * polylog(2.5, -z)

    # Initial guess for the fugacity z.
    # The classical limit corresponds to z -> 0, where the ratio is 1.
    # Since the ratio is 0.75, we expect z to be a positive number.
    # A reasonable initial guess is 1.0.
    initial_guess = 1.0

    # Use fsolve to find the root of the equation.
    z_solution = fsolve(equation_to_solve, initial_guess)

    # The result is an array, so we take the first element.
    z = z_solution[0]

    # Print the final equation with the solved value of z
    f_3_2_val = -polylog(1.5, -z)
    f_5_2_val = -polylog(2.5, -z)
    
    print(f"The problem reduces to solving the equation: f_3/2(z) = 0.75 * f_5/2(z)")
    print(f"Using the numerically found value of z = {z:.2g}:")
    print(f"{f_3_2_val:.4f} = 0.75 * {f_5_2_val:.4f}")
    print(f"{f_3_2_val:.4f} = {0.75 * f_5_2_val:.4f}")
    print("\nThe required value of the fugacity is:")
    # Format the final answer to two significant digits.
    print(f"{z:.2g}")

find_fugacity()