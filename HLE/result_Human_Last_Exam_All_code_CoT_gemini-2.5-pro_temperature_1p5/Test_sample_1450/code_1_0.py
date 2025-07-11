import numpy as np
from scipy.special import polylog
from scipy.optimize import root_scalar

def solve_fugacity():
    """
    Solves for the fugacity z of a Fermi gas where its number density is
    75% of a classical ideal gas at the same pressure and temperature.
    This requires solving the equation: f_{3/2}(z) = 0.75 * f_{5/2}(z),
    which is equivalent to polylog(1.5, -z) - 0.75 * polylog(2.5, -z) = 0.
    """
    
    # Define the equation to be solved for z.
    # The fugacity z must be a positive number.
    def equation(z):
        # The equation is based on f_n(z) = -polylog(n, -z)
        return polylog(1.5, -z) - 0.75 * polylog(2.5, -z)

    # Numerically find the root of the equation.
    # A bracket can be found by testing: equation(0.5) is negative, equation(1.0) is positive.
    # The root must lie between 0.5 and 1.0.
    try:
        solution = root_scalar(equation, bracket=[0.1, 5.0])
        z = solution.root
        
        # The numbers in the final equation f_{n1}(z) = C * f_{n2}(z) are:
        n1 = 1.5
        C = 0.75
        n2 = 2.5
        
        print(f"The equation relating the Fermi-Dirac integrals is: f_{{{n1}}}(z) = {C} * f_{{{n2}}}(z)")
        
        # Print the final result formatted to two significant digits.
        print(f"The value of the fugacity z for which this holds is: {z:.2g}")

    except (ImportError, ValueError) as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have scipy installed (`pip install scipy`).")

if __name__ == "__main__":
    solve_fugacity()
