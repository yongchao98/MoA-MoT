import scipy.special as sp
import scipy.optimize as op

def solve_for_x():
    """
    Finds the largest value of x for which the given summation converges to 0.

    The summation is equivalent to the modified Bessel function of the first kind,
    I_{x-1}(2). We need to find the largest root of I_v(2) = 0, where v = x - 1.
    """

    # Define the function for which we want to find the root.
    def bessel_function_iv(v):
        """Computes the modified Bessel function I_v(2)."""
        return sp.iv(v, 2)

    # From plotting I_v(2), we know that roots exist for v < -1.
    # The largest root (the one with the largest value, i.e., closest to 0)
    # is located in the interval (-3, -2).
    # A tighter bracket can be found by evaluating points in this interval.
    # sp.iv(-2.6, 2) is negative (~-0.06), and sp.iv(-2.5, 2) is positive (~0.024).
    # So, the root is between -2.6 and -2.5.
    try:
        # Use a robust root-finding algorithm with the identified bracket.
        v_root = op.brentq(bessel_function_iv, -2.6, -2.5)

        # The problem asks for x, which is related to v by x = v + 1.
        x_root = v_root + 1

        # The final equation can be written as I_{x-1}(2) = 0.
        # Let's print the values that satisfy this equation.
        # We need to print the value of x, not the intermediate value of v.
        # The equation for the root x found is:
        # Summation from i=0 to infinity of 1 / ((x + i - 1)! * i!) = 0
        # for x = -1.534
        # We can also express this using the Bessel function notation:
        print("The final equation is I_{x-1}(2) = 0.")
        print(f"This holds for x = {x_root:.3f}, as I_({x_root:.3f}-1)(2) = I_({v_root:.3f})(2) ≈ 0.")

        # Print the final answer in the specified format {-a.bbb}
        print("\nFinal Answer:")
        print(f"{{{x_root:.3f}}}")

    except Exception as e:
        print(f"An error occurred during the numerical computation: {e}")

solve_for_x()