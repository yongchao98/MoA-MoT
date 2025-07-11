import numpy as np
from scipy.special import gamma

def solve():
    """
    This function calculates the radius of a spherical balloon, y(t), at t=pi/4
    based on a hypothesized solution form that satisfies the initial conditions.
    """

    # The hypothesized solution form is y(t) = C * (cos(t))^(1/2) * (t^4 + 1)^(-1/4).
    # This form satisfies y'(0) = 0.
    
    # Step 1: Determine the constant C from the initial condition y(0).
    # At t=0, our hypothesized y(t) becomes y(0) = C * (cos(0))^(1/2) * (0^4 + 1)^(-1/4) = C * 1 * 1 = C.
    # Therefore, C is equal to the given value for y(0).
    y0_part1 = 128.0
    y0_part2 = 3**(1/6)
    y0_part3 = gamma(2/3)

    C = 1 / (y0_part1 * y0_part2 * y0_part3)

    # Step 2: Calculate y(pi/4) using the determined function y(t).
    t = np.pi / 4

    # Calculate the components of the solution at t=pi/4
    cos_t_term = np.cos(t)**(1/2)
    poly_term = (t**4 + 1)**(-1/4)

    # The full expression for y(pi/4) is C * cos_t_term * poly_term
    result = C * cos_t_term * poly_term

    # As requested, we will output each number in the final equation.
    # The final equation is: y(pi/4) = y(0) * (cos(pi/4))^(1/2) * ((pi/4)^4 + 1)^(-1/4)
    print("Step 1: Determine the constant C from y(0)")
    print(f"y(0) is the inverse of the product of three numbers:")
    print(f"Number 1: {y0_part1}")
    print(f"Number 2 (3^(1/6)): {y0_part2}")
    print(f"Number 3 (Gamma(2/3)): {y0_part3}")
    print(f"So, C = y(0) = {C}")
    print("\nStep 2: Calculate the value of the other terms at t = pi/4")
    print(f"t = pi/4")
    print(f"(cos(pi/4))^(1/2) = {cos_t_term}")
    print(f"((pi/4)^4 + 1)^(-1/4) = {poly_term}")
    print("\nStep 3: Combine all parts to find the final answer")
    print(f"y(pi/4) = {C} * {cos_t_term} * {poly_term}")
    print(f"Final calculated value for y(pi/4) = {result}")

    return result

# Execute the function to find the answer.
final_answer = solve()
print(f"<<<{final_answer}>>>")