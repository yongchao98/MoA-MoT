import math

def solve_trajectory():
    """
    This function calculates the particle's position x(t) at t = 2*sqrt(3).
    """

    # Step 1: Define the constants and initial values from the problem statement.
    # The initial position x(0) is given by a complex formula.
    # x(0) = 3 + 6^(1/3) * (3 - sqrt(3))^(1/3) + 6^(1/3) * (3 + sqrt(3))^(1/3)
    # This can be rewritten as:
    # x(0) = 3 + (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)
    sqrt3 = math.sqrt(3)
    term_under_root1 = 18 - 6 * sqrt3
    term_under_root2 = 18 + 6 * sqrt3

    x0 = 3 + term_under_root1**(1/3) + term_under_root2**(1/3)

    # The time t is given as 2*sqrt(3).
    t = 2 * sqrt3

    # Step 2: Calculate the change in position.
    # The trajectory is x(t) = x(0) + t^2 / 2.
    position_change = t**2 / 2

    # Step 3: Calculate the final position x(t).
    xt = x0 + position_change

    # Step 4: Print the final equation with all its numerical components.
    # The problem asks to output each number in the final equation.
    print(f"The final equation for the trajectory is: x(t) = x(0) + t^2 / 2")
    print(f"The value of the initial position x(0) is: {x0}")
    print(f"The value of the position change t^2 / 2 is: {position_change}")
    print(f"The final position x({t}) is: {xt}")
    print("\nThe final equation with numbers is:")
    print(f"{xt} = {x0} + {position_change}")


solve_trajectory()