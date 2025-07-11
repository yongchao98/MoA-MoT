import numpy as np
from scipy.special import polylog
from scipy.optimize import brentq

def solve_fugacity_problem():
    """
    Solves for the fugacity z of a Fermi gas under the specified conditions.
    """

    # --- Explanation ---
    print("For an ideal Fermi gas, the number density (nu) and pressure (P) are related to the fugacity (z) by:")
    print("  nu_fermi = C * f_{3/2}(z)")
    print("  P / (k_B * T) = C * f_{5/2}(z)")
    print("where C is a constant and f_n(z) is the Fermi-Dirac integral.\n")

    print("For a classical ideal gas under the same pressure P and temperature T, the number density is:")
    print("  nu_classical = P / (k_B * T)")
    print("Using the pressure equation from the Fermi gas, we can write:")
    print("  nu_classical = C * f_{5/2}(z)\n")

    print("The problem states that nu_fermi = 0.75 * nu_classical.")
    print("Substituting the expressions gives the equation to solve for z:")
    print("  C * f_{3/2}(z) = 0.75 * (C * f_{5/2}(z))")
    print("This simplifies to:")
    print("  f_{3/2}(z) = 0.75 * f_{5/2}(z)\n")

    # --- Numerical Solution ---

    # The Fermi-Dirac integral f_n(z) is related to the polylogarithm function Li_n(z) by:
    # f_n(z) = -Li_n(-z)
    # The `scipy.special.polylog` function calculates Li_n(z).
    def equation_to_solve(z):
        """
        Defines the equation f_{3/2}(z) - 0.75 * f_{5/2}(z) = 0 for the root finder.
        """
        f32 = -polylog(1.5, -z)
        f52 = -polylog(2.5, -z)
        return f32 - 0.75 * f52

    # Find the root numerically in a reasonable interval (e.g., 0.1 to 5.0)
    # The root must be positive.
    try:
        # We search for a root in the interval [0.1, 5.0].
        # The signs of the function at the endpoints must differ for brentq.
        # equation_to_solve(0.1) is positive.
        # equation_to_solve(5.0) is negative.
        z_solution = brentq(equation_to_solve, 0.1, 5.0)

        print(f"Solving this equation numerically gives the fugacity z.")

        # --- Output Final Equation with Numbers ---
        f32_val = -polylog(1.5, -z_solution)
        f52_val = -polylog(2.5, -z_solution)

        print("\nChecking the final equation with the calculated fugacity:")
        print(f"  f_{{3/2}}({z_solution:.3f}) = {f32_val:.4f}")
        print(f"  0.75 * f_{{5/2}}({z_solution:.3f}) = 0.75 * {f52_val:.4f} = {0.75 * f52_val:.4f}")
        print("The two sides are equal, confirming our solution.\n")

        # --- Final Answer ---
        print("The value of the fugacity, rounded to two significant digits, is:")
        print(f"{z_solution:.2g}")

    except (ValueError, RuntimeError) as e:
        print(f"An error occurred during root finding: {e}")
        print("Could not find a solution in the specified interval.")

solve_fugacity_problem()