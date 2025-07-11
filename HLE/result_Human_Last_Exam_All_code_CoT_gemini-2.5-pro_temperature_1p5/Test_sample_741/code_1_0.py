import numpy as np
from scipy.special import iv
from scipy.optimize import root_scalar

def find_largest_x():
    """
    Finds the largest value of x for which the given summation converges to 0.
    The summation is equivalent to the modified Bessel function I_{x-1}(2).
    """

    # We need to solve I_{v}(2) = 0, where v = x - 1.
    # We define the function whose root we want to find.
    def bessel_function_for_solver(v):
        """Modified Bessel function I_v(2) for the root finder."""
        return iv(v, 2.0)

    # From mathematical tables and properties of Bessel functions, the largest root
    # (closest to zero) is located between -2.4 and -2.3.
    # We use the root_scalar function to find this root numerically.
    try:
        solution = root_scalar(bessel_function_for_solver, bracket=[-2.4, -2.3])
        v_max = solution.root

        # The corresponding x value is given by x = v + 1.
        x_max = v_max + 1

        # Print the final answer in the requested format {-a.bbb}.
        print(f"{{{x_max:.3f}}}")

    except ImportError:
        print("Error: The 'scipy' library is required to run this code.")
        print("Please install it via 'pip install scipy'.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    find_largest_x()
