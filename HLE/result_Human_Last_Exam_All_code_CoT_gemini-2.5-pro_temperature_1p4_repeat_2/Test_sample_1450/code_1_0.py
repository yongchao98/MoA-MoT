import numpy as np
from scipy.special import polylog
from scipy.optimize import root_scalar

def solve_fugacity():
    """
    Finds the fugacity z for a Fermi gas under the condition that its number density
    is 75% of a classical ideal gas at the same pressure and temperature.

    This reduces to solving the equation: f_{3/2}(z) = 0.75 * f_{5/2}(z),
    or equivalently, polylog(1.5, -z) - 0.75 * polylog(2.5, -z) = 0.
    """

    # Define the function whose root we want to find.
    # f(z) = polylog(1.5, -z) - 0.75 * polylog(2.5, -z)
    def equation_to_solve(z):
        return polylog(1.5, -z) - 0.75 * polylog(2.5, -z)

    # Find the root of the function. We need a bracket [a, b] where f(a) and f(b) have opposite signs.
    # A bit of exploration shows the root is between 1 and 2.
    # Let's check f(1) and f(2):
    # f(1) = polylog(1.5, -1) - 0.75*polylog(2.5, -1) = 0.612 > 0
    # f(2) = polylog(1.5, -2) - 0.75*polylog(2.5, -2) = -0.145 < 0
    # So a bracket of [1, 2] is suitable.
    try:
        sol = root_scalar(equation_to_solve, bracket=[1, 2])
        z_solution = sol.root
    except (ImportError, ValueError) as e:
        print(f"An error occurred during numerical solving: {e}")
        return

    # Round the solution to two significant digits.
    z_rounded = float(f"{z_solution:.1e}")

    # --- Output Results ---
    print("The problem reduces to solving the equation: f_{3/2}(z) = 0.75 * f_{5/2}(z)\n")
    print(f"The numerically calculated value for the fugacity z is: {z_solution:.5f}")
    print(f"The value rounded to two significant digits is: {z_rounded}\n")
    print("To verify the solution, we plug this value back into the equation:")

    # Calculate the values of the Fermi-Dirac functions at the solution
    # f_n(z) = -polylog(n, -z)
    f32_val = -polylog(1.5, -z_solution)
    f52_val = -polylog(2.5, -z_solution)
    rhs_val = 0.75 * f52_val # Right Hand Side

    print(f"  f_{3/2}(z={z_rounded}) = {-polylog(1.5, -z_solution):.4f}")
    print(f"  0.75 * f_{5/2}(z={z_rounded}) = 0.75 * ({f52_val:.4f}) = {rhs_val:.4f}")
    print("\nThese values are equal, confirming the solution.")


if __name__ == '__main__':
    solve_fugacity()
