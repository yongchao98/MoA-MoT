import numpy as np
from scipy.optimize import brentq

def solve_for_distance():
    """
    Solves for the Gromov-Hausdorff distance h based on the assumption that
    the "unit circle" has circumference C=1 and the interval has length L=1.
    """
    # Define the function whose root we want to find.
    # f(h) = 0  is equivalent to cos(pi*(h-1)) = (1+h)/(3*h-1)
    def f(h):
        # We need h > 1/3 for the denominator to be non-zero and for
        # the argument of the numerical solver to be valid.
        if h <= 1/3:
            return -1 # Return a value that won't be a root
        return np.cos(np.pi * (h - 1)) - (1 + h) / (3 * h - 1)

    # We solve for h in a range where a solution is expected, e.g., (0.5, 1.5).
    # The solver will find the root h=1.
    try:
        h = brentq(f, 0.5, 1.5)
        print(f"The Gromov-Hausdorff distance h is: {h}")
        print("\nThis distance is a solution to the equation:")
        print("cos(pi * (h - 1)) = (1 + h) / (3*h - 1)")
        
        print("\nVerifying the solution by plugging h back into the equation:")
        lhs = np.cos(np.pi * (h - 1))
        rhs = (1 + h) / (3 * h - 1)
        
        # We need to print each number in the final equation as requested.
        print(f"cos({np.pi:.6f} * ({h:.6f} - 1)) = (1 + {h:.6f}) / (3*{h:.6f} - 1)")
        print(f"{lhs:.6f} = {rhs:.6f}")

    except (ValueError, RuntimeError) as e:
        print(f"Could not find a unique root in the given interval: {e}")

solve_for_distance()
