import math

def solve_star_angle_ratio():
    """
    Calculates the ratio (1 - cos(theta_14)) / (1 - cos(theta_34)) based on the
    principles of special relativity and the given observational data.
    """
    # The problem asks for the ratio R = (1 - cos(theta_14)) / (1 - cos(theta_34)).
    # Based on the derivation explained above, this ratio is equal to
    # (1 - cos(theta_12)) / (1 - cos(theta_13)).

    # Step 1: Define the angles given for the second frame of reference.
    # Angle between S1 and S2 is a right angle (pi/2 radians or 90 degrees).
    theta_12 = math.pi / 2
    cos_theta_12 = math.cos(theta_12)

    # Angle between S1 and S3 (and S2 and S3) is 3*pi/4 radians or 135 degrees.
    theta_13 = 3 * math.pi / 4
    cos_theta_13 = math.cos(theta_13)

    # Step 2: Calculate the numerator and denominator of the final expression for the ratio.
    numerator = 1 - cos_theta_12
    denominator = 1 - cos_theta_13

    # Step 3: Calculate the final ratio.
    final_ratio = numerator / denominator

    # Step 4: Print the results, showing the numbers in the final equation as requested.
    print("The problem is to find the value of the ratio (1 - cos(theta_14)) / (1 - cos(theta_34)).")
    print("Based on the principles of special relativity, this ratio simplifies to:")
    print("R = (1 - cos(theta_12)) / (1 - cos(theta_13))")
    print("\nUsing the given values:")
    print(f"cos(theta_12) = cos(pi/2) = {cos_theta_12:.4f}")
    print(f"cos(theta_13) = cos(3*pi/4) = {cos_theta_13:.4f}")
    print("\nThe final equation with numbers is:")
    # The f-string formatting {numerator} and {denominator} will insert the exact float values.
    # The expression shows the direct calculation.
    print(f"R = (1 - {cos_theta_12:.1f}) / (1 - ({cos_theta_13:.4f})) = {numerator} / {denominator:.4f}")
    
    # The exact value is 2 - sqrt(2)
    exact_value_str = "2 - sqrt(2)"
    print(f"\nThis simplifies to the exact value of {exact_value_str}.")
    print(f"The final numerical value is: {final_ratio}")


solve_star_angle_ratio()
<<<2 - 1.414213562373095>>>