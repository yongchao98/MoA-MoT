import scipy.special
import scipy.optimize

def solve_for_x():
    """
    Finds the largest x for which the summation converges to 0.
    The summation is equivalent to the modified Bessel function I_{x-1}(2).
    """

    # We need to solve I_nu(2) = 0 for nu = x - 1.
    # The function to find the root of.
    bessel_func = lambda nu: scipy.special.iv(nu, 2)

    # Based on analysis and testing values, the largest root for nu is in (-1.6, -1.5).
    # I_nu(2) is positive for nu > -1.
    # I(-1.5, 2) approx 0.046 > 0
    # I(-1.6, 2) approx -0.093 < 0
    # The bracket for the root-finding algorithm.
    a = -1.6
    b = -1.5

    # Use the Brent method to find the largest root for nu.
    try:
        largest_nu_root = scipy.optimize.brentq(bessel_func, a, b)
    except (ImportError, ValueError) as e:
        print(f"An error occurred: {e}")
        print("Please ensure scipy is installed ('pip install scipy').")
        return

    # Calculate x from the relation x = nu + 1.
    x_value = largest_nu_root + 1

    # Output the numbers involved in the final calculation.
    print(f"The equation to solve is I_nu(2) = 0, where nu = x - 1.")
    print(f"The largest root found for nu is: {largest_nu_root:.5f}")
    print(f"The final equation for x is: x = {largest_nu_root:.5f} + 1")
    
    # Print the final result in the requested format.
    print(f"\nThe largest value of x is: {x_value:.3f}")

solve_for_x()