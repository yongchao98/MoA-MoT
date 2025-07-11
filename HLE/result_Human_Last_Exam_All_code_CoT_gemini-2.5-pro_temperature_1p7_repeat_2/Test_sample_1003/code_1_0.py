import math

def solve_star_angle_problem():
    """
    Calculates the value of (1 - cos(theta_14)) / (1 - cos(theta_34))
    based on the principles of special relativity.
    """

    # As derived in the plan, the desired ratio is equal to an energy ratio E'_3 / E'_1,
    # which can be determined from the angles given in the second frame of reference.
    # The final expression is 1 / (1 - cos(theta_13)), where theta_13 is 3*pi/4.
    # The angle theta_12 = pi/2 is used to establish the relationship but doesn't
    # appear in the final simplified expression for the energy ratio E'_3 / E'_1.
    # To be fully explicit: E'_3 / E'_1 = (1 - cos(theta_12)) / (1 - cos(theta_13))

    # Angle between S1 and S2 is pi/2.
    theta_12 = math.pi / 2
    cos_theta_12 = math.cos(theta_12)

    # Angle between S1 and S3 is 3*pi/4.
    theta_13 = 3 * math.pi / 4
    cos_theta_13 = math.cos(theta_13)

    # The numerator of the ratio E'_3 / E'_1
    numerator_val = 1 - cos_theta_12

    # The denominator of the ratio E'_3 / E'_1
    denominator_val = 1 - cos_theta_13

    # Calculate the final result
    result = numerator_val / denominator_val

    # Print the final equation with the calculated numbers
    print("The value is calculated using the ratio of energies in the second frame, E'_3 / E'_1.")
    print("This ratio is equal to (1 - cos(theta'_12)) / (1 - cos(theta'_13)).")
    print("\nHere is the calculation with the numerical values:")
    print(f"({numerator_val:.4f}) / ({denominator_val:.4f}) = {result:.8f}")

    print("\nBreaking down the denominator:")
    print(f"1 - cos(3*pi/4) = 1 - ({cos_theta_13:.8f}) = {denominator_val:.8f}")

    # For clarity, also show the exact symbolic result.
    print("\nThe exact symbolic result is 1 / (1 - (-sqrt(2)/2)) = 2 / (2 + sqrt(2)) = 2 - sqrt(2).")


solve_star_angle_problem()
print(f"\n<<<{2 - math.sqrt(2)}>>>")