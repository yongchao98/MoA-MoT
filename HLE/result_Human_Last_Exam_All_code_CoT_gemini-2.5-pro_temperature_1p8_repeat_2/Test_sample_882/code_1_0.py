import math

def solve_for_equilibrium():
    """
    This function implements the plan to solve the problem.
    It numerically solves the equilibrium equation and computes the final required value.
    """
    # The derived polynomial equation for x = 1 - p/N is x^24 - 8x + 7 = 0.
    # The numbers (coefficients) in this final equation are 1 (for x^24), -8 (for x^1), and 7 (for x^0).
    print("The final polynomial equation is of the form c24*x^24 + c1*x^1 + c0 = 0, with numbers:")
    print("c24: 1")
    print("c1: -8")
    print("c0: 7")
    print("-" * 20)

    # We solve for the non-trivial root x using fixed-point iteration on x = (x^24 + 7) / 8.
    # An initial guess of 0.9 is appropriate as the root is known to be between 7/8=0.875 and 1.
    x = 0.9

    # Iterate a sufficient number of times to achieve high precision.
    for _ in range(100):
        x_new = (x**24 + 7) / 8
        if abs(x_new - x) < 1e-15:  # Check for convergence
            x = x_new
            break
        x = x_new

    # We have found the root x for the mixed-strategy equilibrium.
    # The quantity to find is floor(10000 * (1-p)).
    # As established in the plan, 1-p = x^24.
    one_minus_p = x**24

    # Calculate the final value.
    value_to_floor = 10000 * one_minus_p

    # The result is the integer part of this value.
    result = math.floor(value_to_floor)
    
    print(f"The root of the equation is x = {x:.8f}")
    print(f"This corresponds to a probability p = {8*(1-x):.8f}")
    print(f"The value of 1-p is equal to x^24 = {one_minus_p:.8f}")
    print(f"The final expression 10000*(1-p) is: {value_to_floor:.4f}")
    print(f"The floor of this value is: {result}")
    
    # Output the final answer in the required format.
    print(f"\n<<<{result}>>>")

solve_for_equilibrium()