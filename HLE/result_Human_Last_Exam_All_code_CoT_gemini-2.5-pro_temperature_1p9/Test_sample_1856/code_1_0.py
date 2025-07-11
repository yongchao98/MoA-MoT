import math

def calculate_angle():
    """
    Calculates the angle of particle C in the lab frame.
    """
    # Given velocity of particle A
    beta_A = 0.95

    # Calculate the Lorentz factor gamma
    gamma_A = 1.0 / math.sqrt(1 - beta_A**2)

    # The simplified formula for the tangent of the angle in the lab frame
    # for an ultra-relativistic particle C emitted at 45 degrees in the rest frame is:
    # tan(theta_C) = 1 / (gamma_A * (1 + sqrt(2) * beta_A))
    
    sqrt_2 = math.sqrt(2)
    
    # Calculate tan(theta_C)
    tan_theta_C = 1.0 / (gamma_A * (1 + sqrt_2 * beta_A))

    # Calculate the angle theta_C in degrees
    theta_C_rad = math.atan(tan_theta_C)
    theta_C_deg = math.degrees(theta_C_rad)
    
    # As requested, output the numbers in the final equation.
    # We will print the values that go into the final calculation of the angle.
    print(f"The calculation is based on the formula: tan(theta_C) = 1 / (gamma_A * (1 + sqrt(2) * beta_A))")
    print(f"beta_A = {beta_A}")
    print(f"gamma_A = {gamma_A:.4f}")
    print(f"sqrt(2) = {sqrt_2:.4f}")
    print(f"tan(theta_C) = 1 / ({gamma_A:.4f} * (1 + {sqrt_2:.4f} * {beta_A})) = {tan_theta_C:.4f}")
    print(f"theta_C = arctan({tan_theta_C:.4f})")
    print(f"Final Angle (degrees) = {theta_C_deg:.3f}")
    
    # The final answer in the required format
    print(f"<<<{theta_C_deg:.3f}>>>")

calculate_angle()