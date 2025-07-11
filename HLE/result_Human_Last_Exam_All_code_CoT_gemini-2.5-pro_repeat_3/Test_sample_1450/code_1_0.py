import numpy as np
from scipy.special import polylog
from scipy.optimize import fsolve

def solve_fugacity():
    """
    Finds the fugacity z for a Fermi gas under specific conditions.

    The problem reduces to solving the equation:
    f_{3/2}(z) = c * f_{5/2}(z)
    where f_n(z) is the Fermi-Dirac integral and c = 0.75.

    The Fermi-Dirac integral is related to the polylogarithm function by:
    f_n(z) = -polylog(n, -z)

    So the equation to solve is:
    -polylog(3/2, -z) = c * (-polylog(5/2, -z))
    or
    -polylog(3/2, -z) + c * polylog(5/2, -z) = 0
    """
    
    # The constant from the problem statement (75%)
    c = 0.75
    
    # Define the function for which we want to find the root z.
    # The function represents f_{3/2}(z) - c * f_{5/2}(z) = 0
    def equation_to_solve(z):
        return -polylog(1.5, -z) + c * polylog(2.5, -z)

    # An initial guess for the value of z. Based on preliminary analysis,
    # the value should be around 1.
    initial_guess = 1.0

    # Use fsolve to find the root of the equation.
    # fsolve returns an array, so we take the first element.
    z_solution = fsolve(equation_to_solve, initial_guess)[0]

    # Print the final equation being solved with the specific number
    print(f"Solving the equation for fugacity (z): f_3/2(z) = {c} * f_5/2(z)")
    
    # Print the result formatted to two significant digits.
    # The '.2g' format specifier handles this.
    print(f"The required value of the fugacity is: {z_solution:.2g}")

solve_fugacity()