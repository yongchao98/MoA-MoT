import math

def solve_star_angle_problem():
    """
    Calculates the value of (1 - cos(theta_14)) / (1 - cos(theta_34))
    based on the principles of special relativity.

    The problem simplifies to finding the ratio of photon energies E'_3 / E'_1
    in the second reference frame. This ratio is determined by the given
    angles in that frame.
    """

    # The problem reduces to solving for the ratio E'_3 / E'_1.
    # From the derivation, we have the relationship:
    # E'_1 = E'_3 * (1 + 1/sqrt(2))
    # Therefore, the ratio E'_3 / E'_1 is 1 / (1 + 1/sqrt(2)).

    # Calculate the value of 1 + 1/sqrt(2)
    one_over_sqrt2 = 1 / math.sqrt(2)
    denominator = 1 + one_over_sqrt2

    # The ratio is 1 / denominator
    ratio = 1 / denominator

    # The simplified exact expression is 2 - sqrt(2)
    val_sqrt2 = math.sqrt(2)
    exact_value_expr = f"2 - sqrt(2)"
    final_value = 2 - val_sqrt2

    print(f"The problem is to find the value of the expression (1 - cos(theta_14)) / (1 - cos(theta_34)).")
    print(f"This expression simplifies to the ratio of observed photon energies E'_3 / E'_1 in the second reference frame.")
    print(f"Using the given angles, this ratio is found to be 1 / (1 + 1/sqrt(2)).")
    print(f"The simplified exact value of this expression is {exact_value_expr}.")
    print(f"The numbers in the final equation are 2 and sqrt({val_sqrt2**2}) = {val_sqrt2:.4f}.")
    print(f"Final calculated value: {final_value}")

solve_star_angle_problem()