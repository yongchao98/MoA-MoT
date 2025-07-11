import math

def solve_star_angle_ratio():
    """
    Calculates the value of (1 - cos(theta_14)) / (1 - cos(theta_34))
    in the second frame of reference.
    
    The logic is explained in the text above. The final result simplifies to
    the ratio: (1 - cos(theta'_12)) / (1 - cos(theta'_13)).
    """
    
    # In the second frame, the angle between S1 and S2 is a right angle.
    theta_12_prime = math.pi / 2
    
    # In the second frame, the angle between S1 and S3 is 3*pi/4.
    theta_13_prime = 3 * math.pi / 4
    
    # Calculate the values for the numerator and denominator of the ratio.
    # Numerator of the final ratio: 1 - cos(theta'_12)
    numerator_val = 1 - math.cos(theta_12_prime)
    
    # Denominator of the final ratio: 1 - cos(theta'_13)
    denominator_val = 1 - math.cos(theta_13_prime)
    
    # The final result is the ratio of these two values.
    result = numerator_val / denominator_val
    
    print("The final value is determined by the ratio (1 - cos(theta'_12)) / (1 - cos(theta'_13)).")
    print("\n--- Calculation Steps ---")
    print(f"Angle theta'_12 = pi/2, so cos(theta'_12) = {math.cos(theta_12_prime):.1f}")
    print(f"Angle theta'_13 = 3*pi/4, so cos(theta'_13) = {math.cos(theta_13_prime):.7f}")
    print("\nFinal Equation:")
    print(f"Value = (1 - {math.cos(theta_12_prime):.1f}) / (1 - ({math.cos(theta_13_prime):.7f}))")
    print(f"Value = {numerator_val} / {denominator_val:.7f}")
    print(f"Value = {result:.7f}")
    
    # The exact symbolic value is 2 - sqrt(2)
    exact_value = 2 - math.sqrt(2)
    print("\nThe exact symbolic value is 2 - sqrt(2).")
    print(f"Symbolic Value check: 2 - sqrt(2) = {exact_value:.7f}")

solve_star_angle_ratio()
<<<2-sqrt(2)>>>