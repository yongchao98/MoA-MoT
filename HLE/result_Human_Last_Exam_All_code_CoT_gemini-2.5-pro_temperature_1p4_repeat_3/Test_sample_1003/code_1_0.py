import math

def solve_star_angle_ratio():
    """
    Solves the problem based on the Lorentz invariance of the 4-momentum scalar product.
    """
    print("Based on the principles of special relativity, the problem asks for a value which can be simplified to a ratio of photon energies in the second reference frame, R = E'_3 / E'_1.")
    print("This ratio can be calculated from the given angles in the second frame.")
    
    # The ratio R = E'_3 / E'_1 is equivalent to (E'_1 * E'_3) / (E'_1^2).
    # From the problem, we derive:
    # E'_1^2 = (1 - (-1/3)) / (1 - cos(theta'_12))
    # E'_1 * E'_3 = (1 - (-1/3)) / (1 - cos(theta'_13))

    # Given angles in the second frame:
    # Angle between S1 and S2 is a right angle.
    cos_theta_12_prime = 0
    
    # Angle between S1 and S3 is 3*pi/4.
    cos_theta_13_prime = -math.sqrt(2) / 2
    
    # The invariant starting value from frame 1 is (1 - (-1/3)) = 4/3.
    invariant_value = 4/3
    
    # Calculate E'_1^2 and E'_1 * E'_3
    E1_squared = invariant_value / (1 - cos_theta_12_prime)
    E1_E3 = invariant_value / (1 - cos_theta_13_prime)

    print("\nThe equation for the final ratio R can be written as:")
    print("R = (E'_1 * E'_3) / (E'_1^2)")
    print(f"R = [({invariant_value:.4f}) / (1 - {cos_theta_13_prime:.4f})] / [({invariant_value:.4f}) / (1 - {cos_theta_12_prime})]")
    
    # This simplifies to R = (1 - cos_theta_12_prime) / (1 - cos_theta_13_prime)
    # But a more direct way from the logic is R = 1 / (1 - cos_theta_13_prime)
    # Let's calculate R = (E1_E3) / (E1_squared)
    final_ratio = E1_E3 / E1_squared
    
    numerator = 1
    denominator = 1 - cos_theta_13_prime
    
    print("\nThis simplifies to a final equation:")
    print(f"R = {numerator} / (1 - ({cos_theta_13_prime:.4f}))")
    print(f"R = {numerator} / {denominator:.4f}")

    # The exact value is 2 - sqrt(2)
    final_value = 2 - math.sqrt(2)
    
    print(f"\nThe calculated value of the ratio is: {final_value}")

solve_star_angle_ratio()

<<<0.5857864376269049>>>