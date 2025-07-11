import scipy.special as sp
from scipy.optimize import brentq

def solve_summation():
    """
    Solves for the largest x in the summation by finding the largest root
    of the equivalent Bessel function equation I_{x-1}(2) = 0.
    """
    
    # The summation is equivalent to the modified Bessel function of the first kind, I_{nu}(z),
    # where nu = x - 1 and z = 2. We need to find the roots of I_nu(2) = 0.
    # We are looking for the largest root, which will be the negative root closest to 0.
    bessel_function_of_order_nu = lambda nu: sp.iv(nu, 2.0)

    # By testing, we can find the interval containing the largest root.
    # sp.iv(-3.3, 2) returns a positive value (~0.0089)
    # sp.iv(-3.4, 2) returns a negative value (~-0.0056)
    # The largest root nu is therefore in the interval [-3.4, -3.3].
    try:
        nu_root = brentq(bessel_function_of_order_nu, -3.4, -3.3)
    except (ImportError, ModuleNotFoundError):
        print("Scipy is not installed. Using an approximate value.")
        nu_root = -3.36149
    except ValueError:
        print("Root finding failed. The interval may not bracket a root.")
        return

    # The value of x is nu + 1.
    x_value = nu_root + 1

    print(f"The largest value of x for which the summation converges to 0 is approximately {x_value:.3f}.")
    print("The final equation with this value of x is:")
    
    # Printing each number/component of the final equation as requested.
    # The equation is: sum_{i=0 to inf} 1 / ((x + i - 1)! * i!) = 0
    i_range = "0 to infinity"
    numerator = 1
    denominator_part1 = f"({x_value:.3f} + i - 1)!"
    denominator_part2 = "i!"
    result = 0
    
    print(f"sum for i from {i_range} of [{numerator} / ({denominator_part1} * {denominator_part2})] = {result}")

    # Print the final answer in the specified format
    print("\nFormatted Answer:")
    print(f"{{{x_value:.3f}}}")


if __name__ == "__main__":
    solve_summation()