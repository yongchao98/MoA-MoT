import math

def calculate_lab_frame_angle():
    """
    Calculates the angle of particle C in the lab frame.
    
    This function follows the plan:
    1. Define the velocity of particle A.
    2. Calculate the Lorentz factor, gamma_A.
    3. Use the derived formula for the tangent of the lab-frame angle.
    4. Compute the angle in degrees and print the results.
    """
    
    # Given velocity of particle A
    beta_A = 0.95
    
    # 1. Calculate the Lorentz factor, gamma_A.
    # gamma = 1 / sqrt(1 - beta^2)
    gamma_A = 1.0 / math.sqrt(1.0 - beta_A**2)

    # 2. Derive the formula for the angle.
    # The angle of C in A's rest frame, theta*, is given by tan(theta*) = p_x*/p_z*.
    # Since p_x* = P_C*/sqrt(2) and p_z* = P_C*/sqrt(2), tan(theta*) = 1, so theta* = 45 degrees.
    # The Lorentz transformation for angles is tan(theta) = sin(theta*)/(gamma * (cos(theta*) + beta)).
    # Substituting theta* = 45 degrees (sin(45)=cos(45)=1/sqrt(2)):
    # tan(theta) = (1/sqrt(2)) / (gamma_A * (1/sqrt(2) + beta_A))
    # This simplifies to: tan(theta) = 1 / (gamma_A * (1 + sqrt(2)*beta_A))
    
    # 3. Calculate the tangent of the lab-frame angle, tan(theta).
    numerator = 1.0
    denominator = gamma_A * (1.0 + math.sqrt(2) * beta_A)
    tan_theta = numerator / denominator

    # 4. Calculate the angle in degrees.
    theta_rad = math.atan(tan_theta)
    theta_deg = math.degrees(theta_rad)

    # Output the numbers used in the final equation as requested.
    print("--- Calculation Breakdown ---")
    print(f"Velocity of parent particle A, \u03B2_A = {beta_A}")
    print(f"Lorentz factor, \u03B3_A = 1 / \u221A(1 - {beta_A}\u00b2) = {gamma_A:.4f}")
    
    # Printing the components of the final equation for tan(theta)
    print("\nThe equation for the tangent of the angle \u03B8 is:")
    print("tan(\u03B8) = 1 / (\u03B3_A * (1 + \u221A2 * \u03B2_A))")
    print("\nPlugging in the numbers:")
    print(f"tan(\u03B8) = {numerator} / ({gamma_A:.4f} * (1 + {math.sqrt(2):.4f} * {beta_A}))")
    print(f"tan(\u03B8) = {tan_theta:.6f}")
    
    # Print the final result.
    print("\n--- Final Answer ---")
    print(f"The angle \u03B8 between particle C and particle A in the lab frame is: {theta_deg:.3f} degrees")

# Execute the function
calculate_lab_frame_angle()