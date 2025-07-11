import numpy as np
from scipy.special import polylog
from scipy.optimize import fsolve

def solve_fugacity():
    """
    Calculates the fugacity z for an ideal Fermi gas under the specified conditions.

    The condition nu_fermi = 0.75 * nu_classical leads to the equation:
    f_{3/2}(z) = 0.75 * f_{5/2}(z)
    where f_n(z) are the Fermi-Dirac integrals.

    This is solved numerically.
    """
    
    # Define the numbers in the final equation
    n1, d1 = 3, 2
    n2, d2 = 5, 2
    ratio = 0.75
    
    print(f"The number density of the Fermi gas is {ratio} times that of the classical gas.")
    print("This leads to the following equation for the fugacity z:")
    print(f"f_{n1}/{d1}(z) = {ratio} * f_{n2}/{d2}(z)")
    print("\nSolving this equation numerically...")

    # The Fermi-Dirac integral f_n(z) is related to the polylogarithm Li_n(x) by
    # f_n(z) = -Li_n(-z).
    # Our equation becomes: -Li_{3/2}(-z) = 0.75 * (-Li_{5/2}(-z))
    # Or: 0.75 * Li_{5/2}(-z) - Li_{3/2}(-z) = 0
    def equation_to_solve(z):
        # We use scipy.special.polylog(n, x) which computes Li_n(x).
        # We need to solve for z in: ratio * f_{5/2}(z) - f_{3/2}(z) = 0
        val = ratio * (-polylog(n2/d2, -z)) - (-polylog(n1/d1, -z))
        return val

    # An initial guess can be found by series expansion, which gives z ~ 1.
    initial_guess = 1.0

    # Use a numerical solver to find the root of the equation.
    # fsolve returns an array, so we take the first element.
    z_solution = fsolve(equation_to_solve, initial_guess)[0]

    # Print the result formatted to two significant digits.
    print(f"\nThe value of the fugacity is: z = {z_solution:.2g}")

if __name__ == "__main__":
    solve_fugacity()