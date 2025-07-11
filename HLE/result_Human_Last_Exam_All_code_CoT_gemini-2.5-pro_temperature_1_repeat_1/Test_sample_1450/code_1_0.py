import numpy as np
from scipy.special import polylog
from scipy.optimize import root_scalar

def solve_fugacity():
    """
    This function calculates the fugacity 'z' for an ideal Fermi gas
    where its number density is 75% of a classical ideal gas under the
    same pressure and temperature.
    """
    # The condition nu_fermi = 0.75 * nu_classical leads to the equation:
    # f_{3/2}(z) / f_{5/2}(z) = 0.75
    # where f_n(z) is the Fermi-Dirac integral.
    # The Fermi-Dirac integral is related to the polylogarithm function Li_n(x) by:
    # f_n(z) = -Li_n(-z)
    # So, the equation becomes: -polylog(3/2, -z) / -polylog(5/2, -z) = 0.75
    
    ratio = 0.75
    n1 = 1.5  # Corresponds to 3/2
    n2 = 2.5  # Corresponds to 5/2

    def equation_to_solve(z):
        """
        Defines the equation f(z) = 0 that we need to solve for the fugacity z.
        The function should return 0 when the condition is met.
        """
        # polylog(s, z) from scipy.special computes the polylogarithm.
        val_n1 = polylog(n1, -z)
        val_n2 = polylog(n2, -z)
        
        # We need to handle the case where z is very small, which could lead to
        # division by zero if both polylog values are zero.
        if val_n2 == 0:
            # This happens as z -> 0. The ratio approaches 1.
            return 1.0 - ratio

        return val_n1 / val_n2 - ratio

    # We need to find the root of equation_to_solve(z) = 0.
    # By testing, we can find a bracket for the root.
    # For z=1, the ratio is > 0.75.
    # For z=10, the ratio is < 0.75.
    # So the root lies between 1 and 10.
    try:
        solution = root_scalar(equation_to_solve, bracket=[1, 10], method='brentq')
        fugacity = solution.root

        # Output the results as requested
        print(f"The relationship between the number densities leads to the equation for fugacity (z):")
        print(f"f_{n1}(z) / f_{n2}(z) = {ratio}")
        print(f"\nSolving this equation numerically gives the fugacity.")
        print(f"The value of the fugacity is {fugacity:.2g}.")

    except (ValueError, RuntimeError) as e:
        print(f"An error occurred during the numerical solution: {e}")

# Run the solver and print the output
solve_fugacity()