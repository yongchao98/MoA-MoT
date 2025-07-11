import numpy as np
from scipy.special import iv
from scipy.optimize import root_scalar

def solve_summation():
    """
    Finds the largest x for which the summation converges to 0.
    The summation is equivalent to the modified Bessel function of the first kind, I_{x-1}(2).
    We need to find the largest root of I_v(2) = 0, where v = x - 1.
    """

    # Define the function whose root we want to find.
    # This is the modified Bessel function I_v(2).
    def bessel_func_of_order(v):
        return iv(v, 2)

    # We need to find the largest root of I_v(2) = 0.
    # The roots are all negative. We look for the one closest to 0.
    # By plotting or testing, we can find a bracket for the root.
    # iv(-2.7, 2) is positive and iv(-2.8, 2) is negative.
    bracket = [-2.8, -2.7]

    # Use a numerical solver to find the root v.
    try:
        sol = root_scalar(bessel_func_of_order, bracket=bracket, method='brentq')
        v_root = sol.root
    except (ImportError, ValueError) as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have scipy installed (`pip install scipy`).")
        return

    # Calculate the corresponding x value from the relation v = x - 1.
    x_val = 1 + v_root

    # As requested, output the numbers in the final equation: x = 1 + v
    print(f"The summation converges to the modified Bessel function I_v(2), where v = x-1.")
    print(f"We need to find the largest root of I_v(2) = 0.")
    print(f"The largest root is found to be v = {v_root:.10f}")
    print(f"The final equation is x = 1 + v.")
    print(f"So, x = 1 + ({v_root:.3f}) = {x_val:.3f}")

    # The problem asks to respond in the format {-a.bbb}.
    print(f"\nThe largest x value is {{{x_val:.3f}}}")


if __name__ == "__main__":
    solve_summation()