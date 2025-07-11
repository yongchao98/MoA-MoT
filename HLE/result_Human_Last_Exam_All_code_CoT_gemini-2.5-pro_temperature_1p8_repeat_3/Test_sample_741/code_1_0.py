import numpy as np
from scipy.special import iv
from scipy.optimize import brentq

def solve_for_x():
    """
    Finds the largest value of x for which the summation converges to 0.

    The summation is equivalent to the modified Bessel function I_{x-1}(2).
    We need to find the largest root of the equation I_v(2) = 0, where v = x - 1.
    The largest root for v is known to be in the interval (-1.5, -1.0).
    """
    try:
        # Define the function whose root we want to find
        bessel_func = lambda v: iv(v, 2)

        # Find the largest root v_max for I_v(2) = 0 in the interval (-1.5, -1.0)
        v_max = brentq(bessel_func, -1.5, -1.0)

        # The corresponding x value is v_max + 1
        x_max = v_max + 1

        # Print the result in the required format {-a.bbb}
        print(f"{{{x_max:.3f}}}")

    except ImportError:
        print("Error: The 'scipy' library is required.")
        print("Please install it using: pip install scipy numpy")
    except ValueError as e:
        print(f"An error occurred during calculation: {e}")

if __name__ == '__main__':
    solve_for_x()