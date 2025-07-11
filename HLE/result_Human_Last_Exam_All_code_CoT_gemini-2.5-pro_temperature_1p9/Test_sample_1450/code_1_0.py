import scipy.optimize
import scipy.special

def solve_for_fugacity():
    """
    Finds the fugacity z for a Fermi gas where its number density is 75%
    of a classical gas under the same pressure and temperature.

    This involves numerically solving the equation:
    f_{3/2}(z) = 0.75 * f_{5/2}(z)
    where f_n(z) are the Fermi-Dirac integrals.
    """

    # Define the equation whose root (z) we need to find.
    # The equation is f_{3/2}(z) - 0.75 * f_{5/2}(z) = 0.
    # Using the relationship f_n(z) = -polylog(n, -z), we rewrite the equation as:
    # -polylog(3/2, -z) - 0.75 * (-polylog(5/2, -z)) = 0
    # Or: polylog(3/2, -z) - 0.75 * polylog(5/2, -z) = 0
    def equation(z):
        order_1 = 1.5  # This is 3/2
        order_2 = 2.5  # This is 5/2
        factor = 0.75
        # scipy.special.polylog(n, x) is the polylogarithm function Li_n(x).
        return -scipy.special.polylog(order_1, -z) - factor * (-scipy.special.polylog(order_2, -z))

    # The equation to be solved has the following numbers:
    val_1 = 3/2
    val_2 = 0.75
    val_3 = 5/2
    print(f"Solving the equation f_{val_1}(z) = {val_2} * f_{val_3}(z) numerically for z.")

    # Find the root of the equation. Based on physical reasoning and preliminary
    # analysis, the fugacity z is expected to be a positive value around 1.
    # We use a bracket [0.1, 2.0] for the root-finding algorithm.
    try:
        solution = scipy.optimize.root_scalar(equation, bracket=[0.1, 2.0], method='brentq')
        fugacity = solution.root
        # Format the result to two significant digits as requested.
        print(f"The calculated value of the fugacity is: {fugacity:.2g}")
    except (ImportError, ValueError) as e:
        print("An error occurred. Please ensure scipy is installed.")
        print(e)


if __name__ == "__main__":
    solve_for_fugacity()