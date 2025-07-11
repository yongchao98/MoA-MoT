import scipy.special as sp
import scipy.optimize as opt
import numpy as np

def find_largest_x():
    """
    This function finds the largest value of x for which the summation
    S(x) = sum_{i=0 to inf} 1 / ((x + i - 1)! * i!)
    converges to 0.

    The summation S(x) is equivalent to the modified Bessel function of the first kind,
    I_{x-1}(2).

    Thus, we need to find the largest root of the equation I_{nu}(2) = 0, where nu = x - 1.
    """

    # Define the function whose root we want to find.
    # nu represents the order of the Bessel function.
    def bessel_function_at_2(nu):
        return sp.iv(nu, 2)

    # We know the largest root for nu is negative.
    # By testing values, we can find a bracket for the root:
    # I_{-3.0}(2) is positive and I_{-3.5}(2) is negative.
    # Therefore, the largest root lies in the interval [-3.5, -3.0].
    bracket_for_nu = [-3.5, -3.0]

    # Use the Brentq root-finding algorithm to find nu.
    try:
        root_nu = opt.brentq(bessel_function_at_2, bracket_for_nu[0], bracket_for_nu[1])
    except (ValueError, RuntimeError) as e:
        print(f"Root finding failed: {e}")
        return

    # Calculate x from the found root nu, since nu = x - 1.
    root_x = root_nu + 1

    # Print the result in the specified format {-a.bbb}
    print(f"The final equation is: sum_{i=0 to inf} 1 / (({root_x:.3f} + i - 1)! * i!) = 0")
    print(f"The largest value of x for which this holds is approximately: {root_x:.3f}")
    # The final answer format is specific, so we print it exactly as requested.
    print("\nFormatted answer:")
    print(f"{{-{abs(root_x):.3f}}}")

if __name__ == "__main__":
    find_largest_x()
