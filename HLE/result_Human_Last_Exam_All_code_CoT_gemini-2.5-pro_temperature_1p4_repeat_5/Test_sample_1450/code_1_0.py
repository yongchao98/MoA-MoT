import numpy as np
from scipy.special import polylog
from scipy.optimize import root_scalar

def find_fugacity():
    """
    Finds the fugacity z for a Fermi gas where its number density is 75%
    of a classical gas under the same P and T.

    This requires solving the equation:
    f_{3/2}(z) = 0.75 * f_{5/2}(z)
    where f_n(z) is the Fermi-Dirac integral.

    This is equivalent to solving for the root of the function:
    polylog(1.5, -z) - 0.75 * polylog(2.5, -z) = 0
    """

    # The numbers in the final equation are 3/2, 0.75, and 5/2.
    order_f32 = 3.0 / 2.0
    factor = 0.75
    order_f52 = 5.0 / 2.0

    print(f"For an ideal Fermi gas, the number density (nu) and pressure (P) are related to the fugacity (z) by:")
    print(f"nu = C * f_{order_f32}(z)")
    print(f"P/(kT) = C * f_{order_f52}(z)")
    print(f"where C is a constant depending on temperature and particle mass.")
    print("\nFor a classical ideal gas, nu_cl = P/(kT).")
    print(f"The condition nu = {factor} * nu_cl leads to the equation:")
    print(f"f_{order_f32}(z) = {factor} * f_{order_f52}(z)")

    # Define the function whose root we want to find.
    def equation_to_solve(z):
        return polylog(order_f32, -z) - factor * polylog(order_f52, -z)

    # We expect the root to be slightly greater than 1.
    # We can test that f(1) is negative and f(1.5) is positive,
    # so a root exists in the bracket [1.0, 1.5].
    try:
        solution = root_scalar(equation_to_solve, bracket=[1.0, 1.5])
        z_value = solution.root
        print(f"\nSolving this equation numerically gives the fugacity z.")
        print(f"The calculated value of z is: {z_value}")
        print(f"The value rounded to two significant digits is: {z_value:.2g}")
    except (ImportError, ValueError) as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have scipy installed (`pip install scipy`).")

if __name__ == "__main__":
    find_fugacity()
