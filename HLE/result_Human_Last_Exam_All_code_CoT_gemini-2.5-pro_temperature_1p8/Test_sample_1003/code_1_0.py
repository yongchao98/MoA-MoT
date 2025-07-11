import math

def solve_star_angle_ratio():
    """
    Calculates the ratio (1 - cos(theta_14)) / (1 - cos(theta_34))
    based on the principles of special relativity.
    """
    
    # As derived in the plan, the problem simplifies to calculating the ratio
    # of two terms based on the given angles in the second frame of reference.
    # Ratio = (1 - cos(theta'_12)) / (1 - cos(theta'_13))

    # Angles given in the problem for the second frame
    theta_12 = math.pi / 2
    theta_13 = 3 * math.pi / 4

    # Calculate the cosine of these angles
    cos_theta_12 = math.cos(theta_12)
    cos_theta_13 = math.cos(theta_13)

    # Calculate the numerator and denominator of the final expression
    numerator = 1 - cos_theta_12
    denominator = 1 - cos_theta_13

    # Calculate the final result
    result = numerator / denominator
    
    print("The problem is to find the value of (1 - cos(theta_14)) / (1 - cos(theta_34)).")
    print("This ratio simplifies to (1 - cos(theta_12)) / (1 - cos(theta_13)).\n")
    print("Let's calculate the values for the final equation:")
    print(f"cos(theta_12) = cos(pi/2) = {cos_theta_12:.5f}")
    print(f"cos(theta_13) = cos(3*pi/4) = {cos_theta_13:.5f}")
    
    print("\nThe final equation is:")
    # The explicit format string is used to correctly show the double negative.
    print(f"(1 - {cos_theta_12:.5f}) / (1 - ({cos_theta_13:.5f})) = {numerator:.5f} / {denominator:.5f}")

    print(f"\nThe final value is: {result:.5f}")

    # For verification, the exact symbolic result is 2 - sqrt(2)
    exact_result = 2 - math.sqrt(2)
    print(f"The exact value is 2 - sqrt(2), which is approximately {exact_result:.5f}.")

solve_star_angle_ratio()