import scipy.special
import scipy.optimize

def solve_for_x():
    """
    Finds the largest value of 'x' for which the summation
    sum_{i=0 to inf} 1 / ((x + i - 1)! * i!)
    converges to 0.

    This summation is equivalent to the modified Bessel function of the first kind, I_{x-1}(2).
    We need to find the largest root 'v' of the equation I_v(2) = 0, where v = x - 1.
    The roots of I_v(2) with respect to v are all real and negative. The largest root is
    the one closest to zero.
    """
    
    # Define the function for which we are finding the root.
    def bessel_function_of_order_v(v):
        """Computes I_v(2)"""
        return scipy.special.iv(v, 2)

    try:
        # The largest root of I_v(2)=0 is known to be in the interval [-2.8, -2.7].
        # We can confirm this as bessel_function_of_order_v(-2.8) and
        # bessel_function_of_order_v(-2.7) have opposite signs.
        lower_bound = -2.8
        upper_bound = -2.7

        # Use Brent's method to find the root v = x - 1.
        v_root = scipy.optimize.brentq(bessel_function_of_order_v, lower_bound, upper_bound)
        
        # Calculate x from the relation x = v + 1.
        x_root = v_root + 1
        
        # The final equation can be considered x = solution.
        # We print the solution formatted to three decimal places.
        print(f"{x_root:.3f}")

    except ImportError:
        print("Error: This script requires the 'scipy' library.")
        print("Please install it using: pip install scipy")
    except ValueError as e:
        print(f"Error during root finding. The root may not be in the interval [{lower_bound}, {upper_bound}].")
        print(f"Details: {e}")

if __name__ == "__main__":
    solve_for_x()