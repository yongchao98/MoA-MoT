import math

def solve_star_angle_ratio():
    """
    Calculates the ratio (1 - cos(theta_14)) / (1 - cos(theta_34)) based on the
    principles of special relativity and Lorentz invariants.

    The key relation derived is:
    (1 - cos(theta_14)) / (1 - cos(theta_34)) = (1 - cos(theta_12)) / (1 - cos(theta_23))
    """

    # Given values for the second frame of reference
    # Angle between S1 and S2 is a right angle (pi/2 radians)
    cos_theta_12 = 0.0

    # Angle between S2 and S3 is 3*pi/4 radians
    cos_theta_23 = math.cos(3 * math.pi / 4) # This is -1/sqrt(2)

    # Calculate the ratio using the derived formula
    # We handle the case of division by zero, although it won't happen with the given values.
    denominator = (1 - cos_theta_23)
    if denominator == 0:
        print("Error: Division by zero.")
        return

    ratio = (1 - cos_theta_12) / denominator

    # Print the final equation with the numerical values plugged in
    print("The problem reduces to calculating the ratio R = (1 - cos(theta_12)) / (1 - cos(theta_23)).")
    print("Using the given values:")
    print(f"cos(theta_12) = {cos_theta_12}")
    print(f"cos(theta_23) = {cos_theta_23:.10f}")
    print("\nThe final equation with these numbers is:")
    print(f"R = (1 - {cos_theta_12}) / (1 - ({cos_theta_23:.10f}))")
    
    # Print the final result
    print("\nThe calculated value of the ratio is:")
    print(f"R = {ratio:.10f}")
    
    # For verification, the exact value is 2 - sqrt(2)
    exact_value = 2 - math.sqrt(2)
    print(f"\nNote: The exact symbolic value is 2 - sqrt(2), which is approximately {exact_value:.10f}.")


solve_star_angle_ratio()