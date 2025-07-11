import math

def solve_star_angle_problem():
    """
    Calculates the value of (1 - cos(theta_14)) / (1 - cos(theta_34))
    based on the principles of special relativity.
    """
    # The problem reduces to calculating the ratio:
    # (1 - cos(theta'_12)) / (1 - cos(theta'_13))
    # where theta'_12 and theta'_13 are the angles in the second frame.

    # Given angles in the second frame of reference
    theta_12 = math.pi / 2
    theta_13 = 3 * math.pi / 4

    # Calculate the cosines of these angles
    cos_theta_12 = math.cos(theta_12)
    cos_theta_13 = math.cos(theta_13)

    # The numerator of the ratio is 1 - cos(theta'_12)
    numerator = 1 - cos_theta_12
    # The denominator of the ratio is 1 - cos(theta'_13)
    denominator = 1 - cos_theta_13

    # Calculate the final result
    result = numerator / denominator

    # Print the steps of the final calculation
    print("The required ratio (1 - cos(theta_14)) / (1 - cos(theta_34)) is equal to the expression (1 - cos(theta_12)) / (1 - cos(theta_13)).")
    print("\nLet's calculate the values:")
    print(f"Angle theta_12 = pi/2, so cos(theta_12) = {cos_theta_12}")
    print(f"Angle theta_13 = 3*pi/4, so cos(theta_13) = {cos_theta_13:.7f}")
    
    print("\nPlugging these into the final equation:")
    # Using f-string formatting to show the numbers in the equation
    print(f"Result = ({numerator} - ({cos_theta_12})) / ({1.0} - ({cos_theta_13:.7f}))")
    print(f"Result = {numerator} / {denominator:.7f}")
    print(f"Result = {result:.7f}")
    
    # Also print the exact symbolic simplification
    # 1 / (1 + 1/sqrt(2)) = sqrt(2)/(sqrt(2)+1) = sqrt(2)*(sqrt(2)-1) = 2 - sqrt(2)
    exact_value_str = "2 - sqrt(2)"
    print(f"\nThe exact value is {exact_value_str}.")

solve_star_angle_problem()
<<<2 - math.sqrt(2)>>>