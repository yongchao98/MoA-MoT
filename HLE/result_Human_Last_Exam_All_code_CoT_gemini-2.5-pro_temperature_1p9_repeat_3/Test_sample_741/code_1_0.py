import scipy.special as sp
import scipy.optimize as opt

def solve_for_x():
    """
    Finds the largest value of x for which the summation converges to 0.
    The summation is equivalent to the modified Bessel function I_{x-1}(2).
    We need to find the largest root of I_v(2) = 0, where v = x - 1.
    """

    # Define the function for which we want to find the root.
    # The variable is the order 'v' of the Bessel function.
    bessel_function_of_order_v = lambda v: sp.iv(v, 2)

    # The properties of the modified Bessel function I_v(z) tell us that for z > 0,
    # the real roots 'v' are all less than -1.
    # By plotting or evaluating the function, we find the largest root is between -3 and -2.
    bracket = [-3, -2]

    # Use a root-finding algorithm to find the root v in the bracket.
    try:
        solution = opt.root_scalar(bessel_function_of_order_v, bracket=bracket, method='brentq')
        v_root = solution.root
    except ValueError:
        print("Error: The function does not have opposite signs at the bracket endpoints.")
        return

    # The final equation to solve is I_{v}(z) = 0.
    # The numbers in this equation are z=2 and the solution v.
    print(f"The summation is equivalent to the function I_(x-1)(z) with z=2.")
    print(f"We solve the equation I_v(2) = 0 to find the largest root for v.")
    
    # We found the largest root v.
    print(f"The largest root found is v = {v_root:.7f}")
    
    # Calculate x from the relation v = x - 1.
    x_value = v_root + 1
    # The numbers in this part of the equation are v and 1.
    print(f"From the relation x = v + 1, we find x.")
    print(f"x = {v_root:.7f} + 1 = {x_value:.7f}")

    # Print the final answer in the requested format.
    print("\nResult:")
    print(f"The largest x value for which the summation converges to 0 is formatted as {{x.xxx}}:")
    print(f"{{{x_value:.3f}}}")


if __name__ == "__main__":
    solve_for_x()
