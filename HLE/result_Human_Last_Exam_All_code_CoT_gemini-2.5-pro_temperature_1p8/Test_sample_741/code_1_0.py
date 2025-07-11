import numpy as np
from scipy.special import iv
from scipy.optimize import root_scalar

def solve_for_x():
    """
    Finds the largest value of x for which the summation converges to 0.
    The summation is equivalent to the modified Bessel function of the first kind I_{x-1}(2).
    """

    # The problem reduces to solving I_{v}(2) = 0, where v = x - 1.
    # We need to find the largest root v of this equation.
    
    # Define the function whose root we want to find.
    def bessel_function_of_order_v(v):
      """Represents the function I_v(2)"""
      return iv(v, 2.0)

    # Based on the properties of the Bessel function, the largest root
    # is known to lie in the interval [-3, -2].
    try:
        # Find the root v within the bracket [-3, -2]
        solution = root_scalar(bessel_function_of_order_v, bracket=[-3, -2])
        v_root = solution.root

        # Calculate x from the relation x = v + 1
        x_value = v_root + 1
        
        # The final equation relates x to the found root v.
        # We output the numbers used in this final calculation.
        print(f"The largest root for v in I_v(2) = 0 is found to be: {v_root:.3f}")
        print("This root is used to find x from the equation: x = v + 1")
        print("Final calculation:")
        # Showing the equation with the numbers plugged in
        print(f"x = {v_root:.3f} + 1 = {x_value:.3f}")
        
        # Responding with the final value of x in the requested format
        print("\nFinal Answer Format:")
        print(f"{{{x_value:.3f}}}")


    except ValueError:
        print("Error: Could not find a root in the specified interval [-3, -2].")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    solve_for_x()
