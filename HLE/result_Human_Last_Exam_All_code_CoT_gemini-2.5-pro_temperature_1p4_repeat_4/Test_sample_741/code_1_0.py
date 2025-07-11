import numpy as np
from scipy.special import iv
from scipy.optimize import brentq

def solve_for_x():
    """
    Finds the largest value of x for which the given summation converges to 0.

    The summation is equivalent to the modified Bessel function I_{x-1}(2).
    Therefore, we need to find the largest root of the equation I_v(2) = 0, where v = x - 1.
    """

    # The function f(v) = I_v(2) for which we want to find the root.
    z_value = 2.0
    def bessel_function_at_z(v):
        return iv(v, z_value)

    # From mathematical properties, roots for v only exist for v < -1.
    # By testing values, we can find an interval containing the largest root.
    # I_v(2) is positive for v=-2.0 and negative for v=-2.1.
    # So, the largest root (closest to -1) is in the interval [-2.1, -2.0].
    lower_bound = -2.1
    upper_bound = -2.0

    # Use the Brent method to find the root within the determined interval.
    try:
        largest_v_root = brentq(bessel_function_at_z, lower_bound, upper_bound)
    except ValueError:
        print("Error: The function does not have opposite signs at the interval endpoints.")
        return

    # Calculate the corresponding value of x.
    largest_x = largest_v_root + 1

    # Output the details of the equation and its solution.
    print("The summation is equivalent to the modified Bessel function I_{x-1}(2).")
    print("The problem reduces to solving the equation: I_v(z) = 0")
    print(f"The numbers in this equation are: v, z={z_value}, and 0.")
    print(f"The largest root found for v is: {largest_v_root:.3f}")
    print(f"This gives the equation: I_({largest_v_root:.3f})({z_value}) = {bessel_function_at_z(largest_v_root):.3f}")
    
    print("\nThe largest value for x is calculated from v using x = v + 1:")
    print(f"x = {largest_v_root:.3f} + 1 = {largest_x:.3f}")
    
    # This will be used for the final answer block.
    # final_answer = f"{{{largest_x:.3f}}}"
    # print(final_answer)

if __name__ == '__main__':
    solve_for_x()
