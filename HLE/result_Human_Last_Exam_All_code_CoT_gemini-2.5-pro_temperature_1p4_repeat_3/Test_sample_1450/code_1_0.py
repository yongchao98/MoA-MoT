import numpy as np
from scipy.special import polylog
from scipy.optimize import fsolve

def solve_fugacity():
    """
    This function calculates the fugacity z for a Fermi gas under the specified conditions.
    """
    # The problem reduces to solving the equation:
    # f_3/2(z) / f_5/2(z) = 0.75
    # where f_n(z) is the Fermi-Dirac integral of order n.

    # We define a function that represents the equation set to zero,
    # i.e., f_3/2(z) - 0.75 * f_5/2(z) = 0, for a numerical solver.
    # The Fermi-Dirac integral f_n(z) is implemented as -polylog(n, -z) in SciPy.
    def equation_to_solve(z):
        if z <= 0:
            # Fugacity z must be positive
            return np.inf
        
        # The numbers in the equation are n1=1.5, n2=2.5, and the ratio=0.75.
        n1 = 1.5
        n2 = 2.5
        ratio_val = 0.75

        f_n1 = -polylog(n1, -z)
        f_n2 = -polylog(n2, -z)
        
        return f_n1 - ratio_val * f_n2

    # An initial guess is needed for the numerical solver. Asymptotic analysis
    # suggests the value is large, so we start the search at 10.0.
    initial_guess = 10.0
    z_solution = fsolve(equation_to_solve, initial_guess)[0]

    # As requested, we print the numbers in the final equation.
    n1_val = 3/2
    n2_val = 5/2
    ratio_val = 0.75
    
    print(f"The condition ν_fermi = {ratio_val} * ν_classical at the same pressure and temperature leads to the equation:")
    print(f"f_{n1_val}(z) / f_{n2_val}(z) = {ratio_val}")
    print("where z is the fugacity and f_n(z) is the Fermi-Dirac integral.")
    
    print(f"\nSolving this equation numerically gives z ≈ {z_solution:.3f}")
    print(f"The value of the fugacity to two significant digits is: {z_solution:.2g}")

solve_fugacity()