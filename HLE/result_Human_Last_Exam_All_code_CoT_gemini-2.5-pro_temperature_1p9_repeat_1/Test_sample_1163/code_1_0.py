import math

def calculate_angular_distance():
    """
    This function calculates the angular distance between two stars based on the given precession problem.
    """
    # Given constants
    precession_period_years = 26000
    axial_tilt_deg = 23.5
    t_A_years = -3000  # Star A was last on the equator 3000 years ago
    t_B_years = 10000 # Star B will be on the equator in 10000 years

    # --- Step 1: Calculate precession phase angles ---
    # The phase angle phi is the amount of precession in degrees from t=0.
    # phi(t) = (t / T) * 360
    phi_A_deg = (t_A_years / precession_period_years) * 360
    phi_B_deg = (t_B_years / precession_period_years) * 360

    # --- Step 2: Solve for the ecliptic coordinates and alpha_C ---
    # From the problem's symmetry and conditions, we can deduce:
    # lambda_B = 0 degrees
    # lambda_A = phi_A + phi_B
    # alpha_C = phi_B
    lambda_B_deg = 0.0
    lambda_A_deg = phi_A_deg + phi_B_deg
    alpha_C_deg = phi_B_deg
    delta_lambda_deg = lambda_A_deg - lambda_B_deg

    # --- Step 3: Calculate the ecliptic latitude (beta) ---
    # beta is derived from the relation: tan(beta) = -tan(epsilon) * cos(alpha_C)
    axial_tilt_rad = math.radians(axial_tilt_deg)
    alpha_C_rad = math.radians(alpha_C_deg)
    
    tan_beta = -math.tan(axial_tilt_rad) * math.cos(alpha_C_rad)
    beta_rad = math.atan(tan_beta)
    beta_deg = math.degrees(beta_rad)

    # --- Step 4: Calculate the angular distance D ---
    # Use the spherical law of cosines for the ecliptic coordinates.
    # cos(D) = sin(beta_A)*sin(beta_B) + cos(beta_A)*cos(beta_B)*cos(lambda_A-lambda_B)
    # Since beta_A = beta_B = beta, this simplifies to:
    # cos(D) = sin(beta)^2 + cos(beta)^2 * cos(delta_lambda)
    
    delta_lambda_rad = math.radians(delta_lambda_deg)
    sin_beta_sq = math.sin(beta_rad)**2
    cos_beta_sq = math.cos(beta_rad)**2
    cos_delta_lambda = math.cos(delta_lambda_rad)

    cos_D = sin_beta_sq + cos_beta_sq * cos_delta_lambda
    
    # The value of cos(D) might be slightly outside [-1, 1] due to precision, so we clamp it.
    if cos_D > 1.0:
        cos_D = 1.0
    elif cos_D < -1.0:
        cos_D = -1.0
        
    D_rad = math.acos(cos_D)
    D_deg = math.degrees(D_rad)

    # --- Step 5: Print the results ---
    print("The angular distance D is calculated using the formula:")
    print("D = arccos(sin(β)² + cos(β)² * cos(Δλ))")
    print("\nCalculated values:")
    print(f"Ecliptic latitude (β): {beta_deg:.4f} degrees")
    print(f"Ecliptic longitude difference (Δλ): {delta_lambda_deg:.4f} degrees")
    
    print("\nPlugging these values into the formula:")
    print(f"D = arccos(sin({beta_deg:.4f})² + cos({beta_deg:.4f})² * cos({delta_lambda_deg:.4f}))")
    print(f"D = arccos({sin_beta_sq:.4f} + {cos_beta_sq:.4f} * {cos_delta_lambda:.4f})")
    print(f"D = arccos({cos_D:.4f})")
    print(f"\nFinal angular distance between the two stars:")
    print(f"{D_deg:.2f} degrees")

if __name__ == '__main__':
    calculate_angular_distance()