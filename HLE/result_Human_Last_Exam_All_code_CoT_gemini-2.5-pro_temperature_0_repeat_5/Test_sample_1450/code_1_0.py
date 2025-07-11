import numpy as np
from scipy.special import polylog
from scipy.optimize import root_scalar

def solve_fugacity():
    """
    Calculates the fugacity z for an ideal Fermi gas under the condition that
    its number density is 75% of a classical ideal gas at the same P and T.
    """
    # The condition simplifies to the equation: f_{3/2}(z) = 0.75 * f_{5/2}(z)
    # where f_n(z) are the Fermi-Dirac integrals.
    # We need to find the root of the function: f_{3/2}(z) - 0.75 * f_{5/2}(z) = 0.
    
    # The Fermi-Dirac integral f_n(z) is related to the polylogarithm Li_n(x) by:
    # f_n(z) = -Li_n(-z)
    # The polylogarithm function is available in scipy.special as polylog(n, x).
    
    def equation_to_solve(z):
        # The equation is f_{1.5}(z) - 0.75 * f_{2.5}(z) = 0
        f_1_5 = -polylog(1.5, -z)
        f_2_5 = -polylog(2.5, -z)
        return f_1_5 - 0.75 * f_2_5

    # We are solving the equation f_{order1}(z) = factor * f_{order2}(z)
    order1 = 1.5
    factor = 0.75
    order2 = 2.5
    
    print(f"Solving the equation for fugacity (z): f_{order1}(z) = {factor} * f_{order2}(z)")

    # Use a numerical root-finding algorithm. We need to provide a bracket
    # [a, b] where the function values at a and b have opposite signs.
    # A quick check shows that a root lies between 0.1 and 2.0.
    try:
        solution = root_scalar(equation_to_solve, bracket=[0.1, 2.0], method='brentq')
        fugacity = solution.root
        
        # Format the result to two significant digits.
        print(f"The calculated value of the fugacity is: {'{:.2g}'.format(fugacity)}")

    except (ImportError, ValueError) as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have scipy installed ('pip install scipy').")

if __name__ == "__main__":
    solve_fugacity()
