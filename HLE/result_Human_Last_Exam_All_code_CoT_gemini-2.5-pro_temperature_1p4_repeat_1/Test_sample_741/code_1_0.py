import scipy.special as sp
import scipy.optimize as opt

def solve_for_x():
    """
    Finds the largest x for which the summation converges to 0.
    The summation is equivalent to the modified Bessel function I(nu, 2),
    where nu = x - 1. We need to find the largest root of I(nu, 2) = 0.
    """
    # Define the function for the root finder, f(nu) = I(nu, 2)
    bessel_function = lambda nu: sp.iv(nu, 2)

    # The largest root for nu is known to be in the interval [-3.5, -3.0]
    # Find the root nu using the Brent-Dekker method
    nu_root = opt.brentq(bessel_function, -3.5, -3.0)

    # The problem asks for x, and we have the relation x = nu + 1
    x_root = nu_root + 1
    
    # Print the equation relating the solved root nu and the final answer x
    print(f"The root nu is {nu_root:.3f}, and x = nu + 1.")
    print(f"The final equation: {nu_root:.3f} + 1 = {x_root:.3f}")
    
    # Print the final answer for x in the specified format {-a.bbb}
    print(f"The largest x value is: {{{x_root:.3f}}}")

if __name__ == "__main__":
    solve_for_x()