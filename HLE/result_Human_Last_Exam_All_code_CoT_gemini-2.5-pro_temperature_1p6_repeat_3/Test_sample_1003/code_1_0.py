import math

def solve_relativistic_angles():
    """
    Calculates the final value based on the derived ratio.
    The required value (1 - cos(theta_14)) / (1 - cos(theta_34)) simplifies to
    (1 - cos(theta_12)) / (1 - cos(theta_13)) based on the physics of relativistic aberration.
    """
    # The value of theta_12 in the second frame is pi/2
    cos_theta_12 = math.cos(math.pi / 2)
    
    # The value of theta_13 in the second frame is 3*pi/4
    cos_theta_13 = math.cos(3 * math.pi / 4)
    
    # The numerator of the ratio we will calculate
    val_12 = 1 - cos_theta_12
    
    # The denominator of the ratio we will calculate
    val_13 = 1 - cos_theta_13
    
    # The final result is the ratio of these two values.
    # Result = val_12 / val_13
    # Result = (1 - 0) / (1 - (-1/sqrt(2))) = 1 / (1 + 1/sqrt(2))
    # This simplifies to 2 - sqrt(2)
    result = 2 - math.sqrt(2)
    
    print("The final equation for the ratio is (1 - cos(pi/2)) / (1 - cos(3*pi/4))")
    print(f"The numbers in this equation are:")
    print(f"In the numerator: 1 - {cos_theta_12:.1f} = {val_12:.1f}")
    print(f"In the denominator: 1 - ({cos_theta_13:.4f}) = {val_13:.4f}")
    
    print("\nThis simplifies to the expression: 2 - sqrt(2)")
    print("The numbers in the final expression are 2 and the square root of 2.")
    print(f"\nThe calculated value is: {result}")

solve_relativistic_angles()