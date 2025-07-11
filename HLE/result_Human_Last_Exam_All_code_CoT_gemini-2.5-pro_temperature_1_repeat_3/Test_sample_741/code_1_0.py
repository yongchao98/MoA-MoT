import numpy as np
from scipy.special import iv
from scipy.optimize import root_scalar

# This script requires the scipy library. You can install it using:
# pip install scipy

def solve_bessel_root():
    """
    The problem reduces to finding the largest x for which the modified Bessel
    function I_{x-1}(2) is zero. Let v = x - 1. We need to find the largest
    (algebraically) root of I_v(2) = 0.
    """
    
    # Define the equation I_v(2) = 0 that we need to solve for v.
    def bessel_equation(v):
        return iv(v, 2.0)

    # The roots v are all negative. We find the largest root (closest to 0)
    # by searching in the interval where it is known to exist.
    # A plot of the function shows the largest root is between -2.5 and -1.5.
    try:
        solution = root_scalar(bessel_equation, bracket=[-2.5, -1.5])
        v_root = solution.root
    except ValueError as e:
        print(f"Root finding failed: {e}")
        print("Could not find a root in the given bracket.")
        return

    # Calculate the corresponding x value from the relation x = v + 1.
    x_value = v_root + 1
    
    # The prompt asks to output each number in the final equation.
    # The final equation is x = v + 1.
    print(f"The largest root of I_v(2) = 0 is v = {v_root:.3f}")
    print(f"The equation to find x is: x = v + 1")
    print(f"Substituting the value of v: x = {v_root:.3f} + 1 = {x_value:.3f}")
    
    # Print the final answer in the requested format {-a.bbb}
    print("\nFinal Answer:")
    print(f"{{{x_value:.3f}}}")

if __name__ == "__main__":
    solve_bessel_root()