import math

def calculate_angle():
    """
    Calculates the angle of particle C in the lab frame based on the provided parameters.
    """
    # Given parameters
    beta_A = 0.95
    # The rest frame angle theta* is 45 degrees, as derived from the four-momentum components
    # (tan(theta*) = p_x* / p_z* = 1)
    theta_star_deg = 45.0
    theta_star_rad = math.radians(theta_star_deg)

    # Calculate the Lorentz factor gamma
    gamma_A = 1 / math.sqrt(1 - beta_A**2)

    # We use the relativistic angle transformation formula:
    # tan(theta_lab) = sin(theta*) / (gamma * (cos(theta*) + beta))

    # Calculate sin and cos of the rest frame angle
    sin_theta_star = math.sin(theta_star_rad)
    cos_theta_star = math.cos(theta_star_rad)

    # Calculate the numerator and denominator of the tan(theta_lab) expression
    numerator = sin_theta_star
    denominator = gamma_A * (cos_theta_star + beta_A)

    # Calculate tan(theta_lab)
    tan_theta_lab = numerator / denominator

    # Calculate theta_lab in radians and then convert to degrees
    theta_lab_rad = math.atan(tan_theta_lab)
    theta_lab_deg = math.degrees(theta_lab_rad)

    # Print the step-by-step calculation
    print(f"The velocity of particle A is beta_A = {beta_A}")
    print(f"The corresponding Lorentz factor is gamma_A = {gamma_A:.3f}")
    print(f"In the rest frame of A, the angle of particle C is theta* = {theta_star_deg} degrees.")
    print("\nThe formula for the angle in the lab frame (theta) is:")
    print("tan(theta) = sin(theta*) / (gamma_A * (cos(theta*) + beta_A))")
    print("\nPlugging in the numbers:")
    # Printing the equation with each number
    print(f"tan(theta) = {sin_theta_star:.3f} / ({gamma_A:.3f} * ({cos_theta_star:.3f} + {beta_A:.2f}))")
    print(f"tan(theta) = {numerator:.3f} / ({gamma_A:.3f} * {cos_theta_star + beta_A:.3f})")
    print(f"tan(theta) = {numerator:.3f} / {denominator:.3f}")
    print(f"tan(theta) = {tan_theta_lab:.5f}")
    print("\nSolving for the angle in degrees:")
    print(f"theta = arctan({tan_theta_lab:.5f})")
    print(f"The final angle is: {theta_lab_deg:.3f} degrees")

    return round(theta_lab_deg, 3)

# Run the calculation and store the final answer
final_angle = calculate_angle()

# Output the final answer in the specified format
print(f"\n<<<{final_angle}>>>")
