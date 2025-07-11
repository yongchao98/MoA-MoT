import math

def calculate_angle_in_lab_frame():
    """
    Calculates the angle of particle C in the lab frame.
    """
    # 1. Define the given parameters
    beta_A = 0.95
    # From the problem description, we found the angle in the rest frame is 45 degrees.
    theta_star_deg = 45.0
    theta_star_rad = math.radians(theta_star_deg)

    # 2. Calculate intermediate values
    # Calculate the Lorentz factor gamma for particle A
    gamma_A = 1.0 / math.sqrt(1 - beta_A**2)
    
    # Get the sine and cosine of the angle in the rest frame
    sin_theta_star = math.sin(theta_star_rad)
    cos_theta_star = math.cos(theta_star_rad)

    # 3. Apply the relativistic angle transformation formula
    # tan(theta_lab) = sin(theta*) / (gamma * (cos(theta*) + beta))
    tan_theta_lab_numerator = sin_theta_star
    tan_theta_lab_denominator = gamma_A * (cos_theta_star + beta_A)
    tan_theta_lab = tan_theta_lab_numerator / tan_theta_lab_denominator

    # Calculate the angle in the lab frame in radians
    theta_lab_rad = math.atan(tan_theta_lab)

    # Convert the angle to degrees
    theta_lab_deg = math.degrees(theta_lab_rad)

    # 4. Print the final result in the required format
    # The final equation string shows the calculation with numerical values.
    final_equation = (
        f"theta = arctan( sin({theta_star_deg} deg) / "
        f"( (1/sqrt(1-{beta_A**2})) * (cos({theta_star_deg} deg) + {beta_A}) ) )"
    )
    
    final_equation_with_values = (
        f"theta = arctan( {sin_theta_star:.4f} / "
        f"( {gamma_A:.4f} * ({cos_theta_star:.4f} + {beta_A}) ) )"
    )
    
    print("The angle of particle C in the lab frame is calculated as follows:")
    print(final_equation_with_values)
    print(f"theta = {theta_lab_deg:.3f} degrees")


# Execute the calculation and print the result
calculate_angle_in_lab_frame()
print("\n<<<7.589>>>")