import numpy as np
from scipy.special import polylog
from scipy.optimize import root_scalar

def solve_fugacity():
    """
    This function calculates the fugacity z for a Fermi gas under the specified conditions.

    Problem breakdown:
    1. For a classical ideal gas, the number density is nu_classical = P / (k*T).
    2. For a spin-1/2 ideal Fermi gas, the pressure P and number density nu are given by:
       P / (k*T) = (g / lambda^3) * f_{5/2}(z)
       nu_fermi = (g / lambda^3) * f_{3/2}(z)
       where z is the fugacity, g is spin degeneracy (g=2), lambda is the thermal wavelength,
       and f_n(z) are the Fermi-Dirac integrals.

    3. Under the same P and T, we can express the classical density in terms of Fermi gas quantities:
       nu_classical = (g / lambda^3) * f_{5/2}(z)

    4. The problem states nu_fermi = 0.75 * nu_classical.
       Substituting the expressions from (2) and (3):
       (g / lambda^3) * f_{3/2}(z) = 0.75 * (g / lambda^3) * f_{5/2}(z)

    5. This simplifies to the core equation relating the Fermi-Dirac integrals:
       f_{3/2}(z) = 0.75 * f_{5/2}(z)

    6. This function solves the equation f_{3/2}(z) - 0.75 * f_{5/2}(z) = 0 for z numerically.
       The Fermi-Dirac integral is related to the polylogarithm function Li_n(z) by f_n(z) = -Li_n(-z).
    """

    # Define the equation whose root we want to find.
    # The equation is f_{3/2}(z) - 0.75 * f_{5/2}(z) = 0
    # In terms of polylog: -polylog(1.5, -z) - 0.75 * (-polylog(2.5, -z)) = 0
    def equation(z):
        # The equation contains the numbers 3/2=1.5, 0.75, and 5/2=2.5
        val = -polylog(1.5, -z) - 0.75 * (-polylog(2.5, -z))
        return val

    # We use a numerical root-finding algorithm.
    # A quick check shows the root is between z=1 and z=2.
    try:
        solution = root_scalar(equation, bracket=[1.0, 2.0], method='brentq')
        fugacity = solution.root

        # Format the result to two significant digits.
        # The '.2g' format specifier is suitable for this.
        fugacity_rounded = f"{fugacity:.2g}"

        print("The final equation to solve for fugacity (z) is:")
        print("f_{3/2}(z) = 0.75 * f_{5/2}(z)")
        print("\nSolving numerically...")
        print(f"The calculated value of fugacity is: {fugacity}")
        print(f"The fugacity rounded to two significant digits is: {fugacity_rounded}")

    except (ImportError, ValueError) as e:
        print(f"An error occurred: {e}")
        print("Please ensure scipy is installed ('pip install scipy').")

if __name__ == "__main__":
    solve_fugacity()