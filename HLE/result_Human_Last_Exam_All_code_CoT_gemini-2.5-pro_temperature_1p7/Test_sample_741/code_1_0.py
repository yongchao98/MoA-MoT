from scipy.special import iv
from scipy.optimize import brentq

def solve_summation():
    """
    Finds the largest x for which the summation converges to 0.
    The summation is equivalent to the modified Bessel function I_{x-1}(2).
    """
    # Define the function whose root we want to find.
    # The equation is I_v(2) = 0, where v = x - 1.
    def bessel_function_of_order_v(v):
        return iv(v, 2.0)

    # Find the largest root for v. The roots are all v < -1.
    # By testing, we find the largest root is in the interval [-2.6, -2.5].
    # iv(-2.6, 2) is positive and iv(-2.5, 2) is negative.
    bracket_a = -2.6
    bracket_b = -2.5

    largest_root_v = brentq(bessel_function_of_order_v, bracket_a, bracket_b)

    # The final equation to get x is x = v + 1.
    one = 1.0
    largest_x = largest_root_v + one

    # Output the numbers in the final equation as requested.
    print(f"The largest root for the order v is: {largest_root_v:.4f}")
    print(f"Using the relation x = v + 1, we get:")
    print(f"x = {largest_root_v:.4f} + {one:.1f} = {largest_x:.4f}")
    
    # Print the final answer in the required format {-a.bbb}
    print("\nThe final answer is:")
    print(f"{{ {largest_x:.3f} }}")

solve_summation()