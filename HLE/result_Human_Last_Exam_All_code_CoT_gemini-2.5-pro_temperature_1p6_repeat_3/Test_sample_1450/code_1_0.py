import scipy.special
import scipy.optimize

def solve_fugacity():
    """
    Calculates the fugacity z for a Fermi gas where its number density
    is 75% of a classical ideal gas under the same P and T.

    This involves numerically solving the equation:
    f_{3/2}(z) = 0.75 * f_{5/2}(z)

    which is equivalent to:
    polylog(1.5, -z) - 0.75 * polylog(2.5, -z) = 0
    """

    # The equation we need to find the root for is derived from the physics principles.
    # We are solving for z in: Li_{1.5}(-z) - 0.75 * Li_{2.5}(-z) = 0
    ratio = 0.75
    order_nu = 1.5
    order_p = 2.5

    def equation_to_solve(z):
        """ The function for which we are finding the root. """
        return scipy.special.polylog(order_nu, -z) - ratio * scipy.special.polylog(order_p, -z)

    # By testing, we can find that the root lies between z=8 and z=9.
    # equation_to_solve(8) is negative and equation_to_solve(9) is positive.
    try:
        # Use a numerical solver to find the root in the interval [1, 20].
        # A tighter bracket like [8, 9] also works.
        fugacity_z = scipy.optimize.brentq(equation_to_solve, 1, 20)

        print(f"The final equation to be solved is: Li_{order_nu}(-z) = {ratio} * Li_{order_p}(-z)")
        # The .2g format provides two significant digits.
        print(f"The value of the fugacity z for which this holds is: {fugacity_z:.2g}")

    except (ValueError, RuntimeError) as e:
        print(f"An error occurred during root finding: {e}")
        print("Could not find a solution in the given interval.")


if __name__ == '__main__':
    solve_fugacity()
