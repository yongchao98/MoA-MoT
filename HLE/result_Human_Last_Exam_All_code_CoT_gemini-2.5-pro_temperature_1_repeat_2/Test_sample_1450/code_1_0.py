import numpy as np
from scipy.special import polylog
from scipy.optimize import fsolve

def solve_fugacity():
    """
    Calculates the fugacity 'z' for an ideal Fermi gas under the condition
    that its number density is 75% of a classical ideal gas at the same
    pressure and temperature.

    This reduces to solving the equation: 4 * f_{3/2}(z) = 3 * f_{5/2}(z),
    where f_s(z) is the Fermi-Dirac integral.
    """

    # The equation to solve is 4 * f_{3/2}(z) - 3 * f_{5/2}(z) = 0.
    # The Fermi-Dirac integral f_s(z) is implemented as -polylog(s, -z).
    # So, we solve 4 * (-polylog(1.5, -z)) - 3 * (-polylog(2.5, -z)) = 0,
    # which simplifies to 3 * polylog(2.5, -z) - 4 * polylog(1.5, -z) = 0.
    def equation_for_z(z):
        return 3 * polylog(2.5, -z) - 4 * polylog(1.5, -z)

    # An initial guess based on series expansion is z ~ 1.1
    initial_guess = 1.1

    # Use fsolve to find the root of the equation.
    fugacity_solution = fsolve(equation_for_z, initial_guess)

    # Extract the numerical value from the solution array.
    z_value = fugacity_solution[0]
    
    # Define the numbers in the final equation for printing.
    coeff1 = 4
    order1 = 1.5
    coeff2 = 3
    order2 = 2.5

    print("The problem requires solving the following equation for the fugacity (z):")
    print(f"{coeff1} * f_{{{order1}}}(z) = {coeff2} * f_{{{order2}}}(z)")
    print("where f_s(z) represents the Fermi-Dirac integral of order s.")
    print("\nSolving this equation numerically gives the fugacity.")
    
    # Print the final result formatted to two significant digits.
    print(f"\nThe value of the fugacity is: {z_value:.2g}")

solve_fugacity()