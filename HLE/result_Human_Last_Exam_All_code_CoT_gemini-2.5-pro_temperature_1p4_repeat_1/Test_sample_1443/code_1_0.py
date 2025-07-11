import math

def display_solution():
    """
    This program presents the final answer to the geometry problem based on the derivation.
    The answer provides the angles between the asymptotes of the specified conic and the line BC.
    """

    print("The angles are given in terms of delta, the angle between line l and line BC.")
    print("Based on the geometric properties of the construction, the asymptotes are parallel and perpendicular to line l.")
    print("-" * 30)

    # The problem asks to output the numbers in the final equations for the angles.
    # Let the required angles be theta_1 and theta_2.
    # Equation 1: theta_1 = 1 * delta + 0
    # Equation 2: theta_2 = 1 * delta + pi/2

    print("The equation for the first angle is:")
    coeff1 = 1
    const1 = 0
    print(f"angle_1 = {coeff1} * delta + {const1}")

    print("\nThe equation for the second angle is:")
    coeff2 = 1
    const2_str = "pi/2"  # Using a string for symbolic representation
    print(f"angle_2 = {coeff2} * delta + {const2_str}")
    print("-" * 30)

display_solution()
