import math

def calculate_decay_angle():
    """
    Calculates the angle between particle C and its parent particle A in the lab frame.
    """
    # Given velocity of particle A
    beta_A = 0.95

    # Calculate the Lorentz factor gamma for particle A
    gamma_A = 1 / math.sqrt(1 - beta_A**2)

    # The formula for the tangent of the angle theta_C in the lab frame is derived from
    # the Lorentz transformation of momentum:
    # tan(theta_C) = 1 / (gamma_A * (1 + sqrt(2) * beta_A))
    
    # Calculate the components of the equation
    sqrt_2 = math.sqrt(2)
    denominator = gamma_A * (1 + sqrt_2 * beta_A)
    tan_theta_C = 1 / denominator

    # Calculate the angle in radians and then convert to degrees
    theta_C_rad = math.atan(tan_theta_C)
    theta_C_deg = math.degrees(theta_C_rad)

    # Output the final equation with all numbers plugged in
    print("The final equation for the tangent of the angle is:")
    print(f"tan(theta_C) = 1 / ({gamma_A:.5f} * (1 + {sqrt_2:.5f} * {beta_A}))")
    print(f"tan(theta_C) = 1 / ({denominator:.5f}) = {tan_theta_C:.5f}")
    
    # Print the final result
    print("\nFinal Answer:")
    print(f"The angle between the direction of C and that of its parent particle A is {theta_C_deg:.3f} degrees.")

calculate_decay_angle()