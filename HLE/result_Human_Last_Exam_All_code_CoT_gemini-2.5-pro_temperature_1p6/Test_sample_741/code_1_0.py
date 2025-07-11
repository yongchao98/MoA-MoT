import scipy.special as sp
from scipy.optimize import brentq

def solve_for_x():
    """
    This function finds the largest value of x for which the summation
    Sum_{i=0 to inf} 1 / ((x + i - 1)! * i!)
    converges to 0.

    The summation is equivalent to the modified Bessel function of the first kind, I_{x-1}(2).
    Therefore, we need to solve the equation I_{x-1}(2) = 0.
    Let v = x - 1. We find the largest root of I_v(2) = 0.
    """

    # Define the function whose root we want to find. f(v) = I_v(2)
    def bessel_function_iv(v):
        return sp.iv(v, 2.0)

    # The roots of I_v(2) = 0 are all negative. From numerical analysis,
    # the largest root (the one closest to 0) lies in the interval [-2.8, -2.7].
    # We use the brentq root-finding algorithm on this interval.
    # Note: bessel_function_iv(-2.8) is negative and bessel_function_iv(-2.7) is positive.
    try:
        # Find the largest root for v.
        v_root = brentq(bessel_function_iv, -2.8, -2.7)
    except ImportError:
        print("Scipy library not found. Please install it using 'pip install scipy'")
        return
    except ValueError as e:
        print(f"Root finding failed. The interval might not bracket a root. Error: {e}")
        return


    # The value of x is v + 1.
    x_value = v_root + 1

    # Print the result in the specified format {-a.bbb}
    print(f"{{{x_value:.3f}}}")

    # Print the final equation with the calculated numbers
    print("\nThe equation for which this value of x is the solution:")
    print(f"Sum_{{i=0 to inf}} 1 / (({x_value:.3f} + i - 1)! * i!) = 0")
    print("Which simplifies to:")
    print(f"Sum_{{i=0 to inf}} 1 / ((i + {v_root:.3f})! * i!) = 0")


if __name__ == '__main__':
    solve_for_x()