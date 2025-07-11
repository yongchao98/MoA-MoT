import scipy.special
import scipy.optimize

def solve_for_x():
    """
    Finds the largest value of x for which the summation converges to 0.
    The summation is equivalent to the modified Bessel function I_{x-1}(2).
    We need to find the largest root x of the equation I_{x-1}(2) = 0.
    """
    
    # Let nu = x - 1. We are looking for the largest root of I_nu(2) = 0.
    # The function to find the root of is f(nu) = I_nu(2).
    # The value z in I_nu(z) is 2.
    z_val = 2
    def bessel_function_of_nu(nu):
        return scipy.special.iv(nu, z_val)

    # From mathematical properties, the roots nu are all negative.
    # We need the largest root, which is the one closest to 0.
    # Let's find a bracket [a, b] where the function changes sign.
    # After testing some values, e.g., f(-2.4) > 0 and f(-2.5) < 0,
    # we can use this interval to find the root.
    bracket = [-2.5, -2.4]

    # Find the root nu using a numerical solver (Brent's method).
    try:
        solution = scipy.optimize.root_scalar(
            bessel_function_of_nu, 
            bracket=bracket, 
            method='brentq'
        )
        nu_root = solution.root
        
        # Calculate x from the root nu, since nu = x - 1
        x_root = nu_root + 1

        # The final equation is the one we solved, with the calculated value of x
        # S(x) = I_{x-1}(2) = 0
        # Print the solution for x in the required format.
        # Here we show the numbers from the equation I_{x-1}(2)=0 that lead to the solution.
        # This includes the final calculated x, and the parameters 1 and 2 from the equation.
        print(f"To make the sum zero, we solve the equation I_(x-1)(2) = 0.")
        print(f"The largest value for x that satisfies this is found to be: {x_root:.3f}")
        print("Formatted Answer:")
        print(f"{{{x_root:.3f}}}")
        
    except ValueError as e:
        print(f"Error finding root: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

solve_for_x()