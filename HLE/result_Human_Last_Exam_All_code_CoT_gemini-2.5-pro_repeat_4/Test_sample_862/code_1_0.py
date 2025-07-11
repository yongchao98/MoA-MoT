import numpy as np
from scipy.optimize import brentq

def main():
    """
    This script solves for the constant C by finding the root of a
    transcendental equation and then using it to compute C.
    """

    # Define the transcendental equation to be solved
    def equation(x):
        return np.tan(x) - x - np.pi / 2

    # Find the root x0 in the interval (0, pi/2).
    # We use a bracket slightly smaller than (0, pi/2) to avoid the
    # singularity of tan(x) at x = pi/2.
    try:
        x0 = brentq(equation, 0, np.pi / 2 - 1e-9)
    except ValueError:
        print("Could not find a root in the given interval.")
        return

    # Calculate the constant C using the solution x0
    cos_x0 = np.cos(x0)
    C = (1 + cos_x0) / (1 - cos_x0)

    print("The problem reduces to solving the equation tan(x) = x + pi/2.")
    print(f"The numerical solution is x ≈ {x0}")
    print("\nThe constant C is given by the formula C = (1 + cos(x)) / (1 - cos(x)).")
    print("Substituting the value of x found above:")
    
    # We print the final equation with the computed numerical values.
    print(f"C = (1 + {cos_x0}) / (1 - {cos_x0})")
    
    # Print the final result
    print(f"\nThe smallest possible constant C is approximately {C}")

if __name__ == "__main__":
    main()
