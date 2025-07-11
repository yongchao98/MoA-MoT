import numpy as np
from scipy.special import polylog
from scipy.optimize import fsolve

def solve_fugacity():
    """
    Calculates the fugacity z for a Fermi gas under the specified conditions.
    """
    # The problem reduces to solving the transcendental equation:
    # g * f_{3/2}(z) = 0.75 * f_{5/2}(z)
    # For spin-1/2 particles, the degeneracy g = 2.
    # The Fermi-Dirac integral f_k(z) is related to the polylogarithm function Li_k(x)
    # by f_k(z) = -Li_k(-z). The scipy function is `polylog(k, x)`.
    # So, the equation is: 2 * [-polylog(1.5, -z)] = 0.75 * [-polylog(2.5, -z)]
    # which simplifies to 2 * polylog(1.5, -z) = 0.75 * polylog(2.5, -z).

    # We need to find the root of the function F(z) = 0, where:
    # F(z) = 2 * polylog(1.5, -z) - 0.75 * polylog(2.5, -z)

    def equation_to_solve(z):
        """The function whose root we want to find."""
        # The fugacity z must be positive. The polylog function in scipy handles
        # negative arguments, which corresponds to -z.
        return 2.0 * polylog(1.5, -z) - 0.75 * polylog(2.5, -z)

    # An initial guess is needed for the numerical solver. A simple series expansion
    # suggests the root is around z=2.
    initial_guess = 2.0
    
    # Use fsolve to find the root.
    z_solution = fsolve(equation_to_solve, initial_guess)
    
    # fsolve returns an array, so we take the first element.
    z_value = z_solution[0]

    # The final equation being solved is g * f_3/2(z) = C * f_5/2(z)
    g = 2.0
    C = 0.75
    
    print(f"The equation to solve for the fugacity (z) is: {g} * f_3/2(z) = {C} * f_5/2(z)")
    print(f"The numerical solution for z is: {z_value}")
    
    # Format the result to two significant digits.
    # We use the 'g' format specifier for significant figures.
    z_rounded = f"{z_value:.2g}"
    print(f"The value of the fugacity rounded to two significant digits is: {z_rounded}")

solve_fugacity()