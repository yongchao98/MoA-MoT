import math

def calculate_angle():
    """
    Calculates the angle of particle C in the lab frame.
    
    The function follows these steps:
    1. Define the velocity of particle A (beta_A).
    2. Calculate the corresponding Lorentz factor (gamma_A).
    3. Use the relativistic angle transformation formula to find tan(theta), where theta is the angle in the lab frame.
       The formula is derived from the Lorentz transformation of momentum:
       tan(theta) = p_T / p_z
                  = (p_x^*) / (gamma_A * (p_z^* + beta_A * E_C^*))
       With p_x^*=P_C^*/sqrt(2), p_z^*=P_C^*/sqrt(2) and E_C^*=P_C^*, this simplifies to:
       tan(theta) = (1/sqrt(2)) / (gamma_A * (beta_A + 1/sqrt(2)))
    4. Calculate the angle theta in degrees and round to three decimal places.
    5. Print all the intermediate values and the final result.
    """
    
    # Given velocity of particle A
    beta_A = 0.95
    
    # Calculate Lorentz factor gamma_A
    gamma_A = 1 / math.sqrt(1 - beta_A**2)
    
    # Terms in the tan(theta) equation
    one_over_sqrt2 = 1 / math.sqrt(2)
    
    # Calculate tan(theta)
    tan_theta_numerator = one_over_sqrt2
    tan_theta_denominator = gamma_A * (beta_A + one_over_sqrt2)
    tan_theta = tan_theta_numerator / tan_theta_denominator
    
    # Calculate the angle in degrees
    theta_rad = math.atan(tan_theta)
    theta_deg = math.degrees(theta_rad)
    
    # Print the equation and results
    print("Finding the angle (theta) of particle C in the lab frame.\n")
    print(f"The equation for the tangent of the angle is:")
    print("tan(theta) = (1/sqrt(2)) / (gamma_A * (beta_A + 1/sqrt(2)))\n")
    print("Here are the numerical values:")
    print(f"beta_A = {beta_A}")
    print(f"gamma_A = 1 / sqrt(1 - {beta_A}^2) = {gamma_A:.4f}")
    print(f"1/sqrt(2) = {one_over_sqrt2:.4f}")
    
    print(f"\nCalculating tan(theta):")
    print(f"tan(theta) = {tan_theta_numerator:.4f} / ({gamma_A:.4f} * ({beta_A} + {one_over_sqrt2:.4f}))")
    print(f"tan(theta) = {tan_theta_numerator:.4f} / {tan_theta_denominator:.4f}")
    print(f"tan(theta) = {tan_theta:.4f}")
    
    print(f"\nCalculating the angle:")
    print(f"theta = atan({tan_theta:.4f})")
    print(f"theta = {theta_deg:.3f} degrees")
    
if __name__ == "__main__":
    calculate_angle()
    # For the final answer format
    beta_A = 0.95
    gamma_A = 1 / math.sqrt(1 - beta_A**2)
    one_over_sqrt2 = 1 / math.sqrt(2)
    tan_theta = one_over_sqrt2 / (gamma_A * (beta_A + one_over_sqrt2))
    theta_deg = math.degrees(math.atan(tan_theta))
    print(f"\n<<<{theta_deg:.3f}>>>")
