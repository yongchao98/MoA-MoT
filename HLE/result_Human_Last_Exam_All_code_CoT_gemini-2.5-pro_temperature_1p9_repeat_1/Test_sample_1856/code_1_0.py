import math

def calculate_lab_frame_angle():
    """
    Calculates the angle of particle C in the lab frame based on the given parameters.
    """
    # --- Step 1: Define given values and initial conditions ---
    
    # Velocity of the parent particle A in the lab frame
    beta_A = 0.95
    
    # In A's rest frame, the momentum of C has components:
    # p_x* = P_C* / sqrt(2)
    # p_z* = P_C* / sqrt(2)
    # The angle in the rest frame, theta_C*, is given by tan(theta_C*) = p_x*/p_z* = 1.
    # Therefore, theta_C* is 45 degrees.
    theta_C_star_deg = 45.0
    theta_C_star_rad = math.radians(theta_C_star_deg)
    
    # --- Step 2: Calculate intermediate physical quantities ---
    
    # Calculate the Lorentz factor gamma for particle A
    gamma_A = 1 / math.sqrt(1 - beta_A**2)
    
    # Pre-calculate sin and cos of the rest-frame angle
    sin_theta_C_star = math.sin(theta_C_star_rad)
    cos_theta_C_star = math.cos(theta_C_star_rad)
    
    # --- Step 3: Apply the relativistic angle transformation formula ---
    
    # The formula is tan(theta_C) = sin(theta_C*) / (gamma_A * (cos(theta_C*) + beta_A)).
    # We use the simplified version as the condition m_C << E_C* implies the
    # speed of particle C in A's rest frame (beta_C*) is approximately 1.
    
    # Calculate the numerator and denominator of the tangent function
    numerator = sin_theta_C_star
    denominator = gamma_A * (cos_theta_C_star + beta_A)
    
    # Calculate the tangent of the lab frame angle theta_C
    tan_theta_C = numerator / denominator
    
    # --- Step 4: Calculate the final angle and print the results ---
    
    # Calculate the angle theta_C in radians, then convert to degrees
    theta_C_rad = math.atan(tan_theta_C)
    theta_C_deg = math.degrees(theta_C_rad)
    
    # Print the explanation and the equation with numbers filled in as requested
    print("The lab frame angle (theta_C) is found from the rest frame angle (theta_C*)")
    print("using the formula: tan(theta_C) = sin(theta_C*) / (gamma_A * (cos(theta_C*) + beta_A))")
    print("\n--- Calculation Steps ---")
    print(f"1. Angle in A's rest frame, theta_C* = {theta_C_star_deg:.1f} degrees.")
    print(f"   - sin({theta_C_star_deg:.1f}°) = {sin_theta_C_star:.4f}")
    print(f"   - cos({theta_C_star_deg:.1f}°) = {cos_theta_C_star:.4f}")
    print(f"2. Lorentz factor for A, gamma_A = 1 / sqrt(1 - {beta_A}^2) = {gamma_A:.4f}")
    
    print("\n3. Substituting these numbers into the formula:")
    final_equation_str = (
        f"tan(theta_C) = {sin_theta_C_star:.4f} / ({gamma_A:.4f} * ({cos_theta_C_star:.4f} + {beta_A}))"
    )
    print(final_equation_str)
    
    print("\n4. Solving for theta_C:")
    print(f"tan(theta_C) = {sin_theta_C_star:.4f} / {denominator:.4f} = {tan_theta_C:.4f}")
    print(f"theta_C = arctan({tan_theta_C:.4f})")
    print(f"The final angle is: {theta_C_deg:.3f} degrees.")


if __name__ == '__main__':
    calculate_lab_frame_angle()