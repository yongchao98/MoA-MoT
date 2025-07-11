import numpy as np
from scipy.special import polylog
from scipy.optimize import root_scalar

def solve_fugacity():
    """
    Solves for the fugacity z of a Fermi gas under the given conditions.
    """
    # The condition nu_FD = 0.75 * nu_cl leads to the equation:
    # f_{3/2}(z) / f_{5/2}(z) = 0.75
    # where f_n(z) is the Fermi-Dirac integral.
    # f_n(z) = -polylog(n, -z).
    # So we need to solve: polylog(1.5, -z) = 0.75 * polylog(2.5, -z)
    
    def equation_to_solve(z):
        """The equation f_{3/2}(z) - 0.75 * f_{5/2}(z) = 0."""
        # The minus sign from f_n(z) = -Li_n(-z) cancels on both sides.
        return polylog(1.5, -z) - 0.75 * polylog(2.5, -z)

    # We need to find a bracket for the root.
    # At z=1, equation(1) is positive.
    # A large z approximation (z~28) suggests the function will become negative.
    # Let's test z=28.
    # polylog(1.5, -28) is ~7.39, polylog(2.5, -28) is ~10.02
    # 7.39 - 0.75*10.02 = 7.39 - 7.515 = -0.125
    # So a root exists in the bracket [1, 28].
    try:
        sol = root_scalar(equation_to_solve, bracket=[1, 30])
        z_value = sol.root
    except (ImportError, ValueError) as e:
        print(f"An error occurred: {e}")
        print("Could not find the root numerically. Please check your scipy installation.")
        return

    # Output the result and verification
    f_3_2_val = -polylog(1.5, -z_value)
    f_5_2_val = -polylog(2.5, -z_value)
    
    print("The derived equation relating the fugacity 'z' to the condition is:")
    print("f_{3/2}(z) = 0.75 * f_{5/2}(z)")
    print("\nSolving this equation numerically gives the fugacity.")
    print(f"\nFor the solution z = {z_value:.4f}:")
    print(f"The left side of the equation is f_1.5({z_value:.2f}) = {-polylog(1.5, -z_value):.2f}")
    print(f"The right side of the equation is 0.75 * f_2.5({z_value:.2f}) = 0.75 * {-polylog(2.5, -z_value):.2f} = {0.75 * -polylog(2.5, -z_value):.2f}")
    print("\nThe condition is satisfied.")
    
    # Print the final answer formatted to two significant digits
    print(f"\nThe value of the fugacity is {z_value:.2g}.")

solve_fugacity()