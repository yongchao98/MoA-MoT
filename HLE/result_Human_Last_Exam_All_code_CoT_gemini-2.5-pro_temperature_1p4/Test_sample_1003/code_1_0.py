import math

def solve_star_angles():
    """
    Calculates the ratio (1 - cos(theta_14)) / (1 - cos(theta_34))
    based on the provided relativistic observations.
    """
    # Step 1: Define the known cosines of angles from the second frame of reference.
    # Angle between S1 and S2 is pi/2.
    cos_theta_12 = 0
    # Angle between S1/S2 and S3 is 3*pi/4.
    cos_theta_13 = -1 / math.sqrt(2)
    cos_theta_23 = -1 / math.sqrt(2)

    # Step 2: Use the relation 1 - cos(theta_ij) = A_i * A_j to set up equations.
    # From theta_12: A1 * A2 = 1 - 0 = 1
    # From theta_13: A1 * A3 = 1 - (-1/sqrt(2)) = 1 + 1/sqrt(2)
    # From theta_23: A2 * A3 = 1 - (-1/sqrt(2)) = 1 + 1/sqrt(2)
    
    # From A1*A3 = A2*A3, we can deduce A1 = A2 (since A3 cannot be zero).
    # Substituting A1 = A2 into A1*A2 = 1 gives A1^2 = 1.
    # The factor A_i must be positive, so A1 = 1.
    A1 = 1.0
    
    # Calculate A3 from A1 * A3 = 1 + 1/sqrt(2)
    A3 = 1 + 1 / math.sqrt(2)
    
    # Step 3: The target quantity is (1 - cos(theta_14)) / (1 - cos(theta_34)).
    # This translates to (A1 * A4) / (A3 * A4), which simplifies to A1 / A3.
    result = A1 / A3

    # Step 4: Print the results, showing the final equation with the calculated numbers.
    print(f"Based on the given angles, we derive the values for the factors A1 and A3.")
    print(f"A1 = {A1}")
    print(f"A3 = 1 + 1/sqrt(2) = {A3}")
    print("\nThe expression to evaluate is (1 - cos(theta_14)) / (1 - cos(theta_34))")
    print("This simplifies to the ratio A1 / A3.")
    print("\nThe final equation is:")
    print(f"{A1} / {A3}")
    print("\nSymbolically, this is 1 / (1 + 1/sqrt(2)), which is equal to 2 - sqrt(2).")
    print(f"\nThe final numerical value is: {result}")

solve_star_angles()