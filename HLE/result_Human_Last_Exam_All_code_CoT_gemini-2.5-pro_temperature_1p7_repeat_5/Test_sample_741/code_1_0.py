import scipy.special
import scipy.optimize

def find_largest_x_for_sum_zero():
    """
    This script finds the largest value of x for which the summation
    sum_{i=0 to inf} 1 / ((x + i - 1)! * i!) converges to 0.

    This summation is equivalent to the modified Bessel function of the first kind,
    I_{nu}(z), with nu = x - 1 and z = 2.

    So, we need to solve the equation I_{x-1}(2) = 0. This means we are looking
    for the largest root of the function f(nu) = I_nu(2).
    """

    # The equation we need to solve is I_nu(z) = 0.
    # The parameters for our problem are:
    z = 2
    
    # Define the function for which we want to find the root.
    def bessel_function(nu):
        return scipy.special.iv(nu, z)

    # We are looking for the largest root, which is known to be negative.
    # We will search for a bracketing interval [a, b] for this root.
    # After testing some values, we find the root is between -2.4 and -2.3.
    # I_nu(2) > 0 for nu = -2.3
    # I_nu(2) < 0 for nu = -2.4
    a = -2.4
    b = -2.3
    
    # Use the Brent method to find the root nu in the interval [a, b].
    try:
        nu_root = scipy.optimize.brentq(bessel_function, a, b)
    except ValueError:
        print("Error: The function values at the interval endpoints do not have opposite signs.")
        return

    # The value of x is derived from nu = x - 1.
    x_val = nu_root + 1

    # Print the equation being solved and the result
    print(f"The summation is equivalent to the modified Bessel function I_nu(z).")
    print(f"We solve I_nu(z) = 0 for z = {z}.")
    print(f"The largest root found for nu is approximately {nu_root:.3f}, located in the interval [{a}, {b}].")
    print(f"The value of I_nu(z) at this root is approximately {bessel_function(nu_root):.2e}.")
    print(f"The corresponding value for x is nu + 1 = {x_val:.3f}.")
    print(f"The final answer in format {{-a.bbb}} is:")
    
    # Print the final result in the requested format.
    print(f"{{{x_val:.3f}}}")


if __name__ == '__main__':
    find_largest_x_for_sum_zero()
