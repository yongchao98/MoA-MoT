import math

def calculate_decay_angle():
    """
    Calculates the angle of a decay product particle in the lab frame.

    This function follows the plan:
    1. Define the velocity of the parent particle A.
    2. Calculate the Lorentz factor gamma.
    3. Use the derived formula for the tangent of the angle of particle C.
    4. Compute the angle in degrees.
    5. Print the breakdown of the calculation and the final result.
    """
    # 1. Given velocity of particle A
    beta_A = 0.95

    # 2. Calculate the Lorentz factor gamma
    gamma_A = 1 / math.sqrt(1 - beta_A**2)

    # 3. Use the derived formula: tan(theta_C) = (1/sqrt(2)) / (gamma_A * (1/sqrt(2) + beta_A))
    # Let's calculate the components of this equation.
    term_1_over_sqrt_2 = 1 / math.sqrt(2)
    
    # Numerator of the final expression for tan(theta_C)
    numerator = term_1_over_sqrt_2
    
    # Denominator of the final expression for tan(theta_C)
    denominator = gamma_A * (term_1_over_sqrt_2 + beta_A)

    # Calculate the tangent of the angle
    tan_theta_C = numerator / denominator
    
    # 4. Compute the angle in degrees
    angle_rad = math.atan(tan_theta_C)
    angle_deg = math.degrees(angle_rad)

    # 5. Print the equation and the values substituted into it
    print("The final equation for the tangent of the angle is:")
    print("tan(theta_C) = (1/sqrt(2)) / (gamma_A * (1/sqrt(2) + beta_A))")
    print("\nSubstituting the numerical values:")
    print(f"beta_A = {beta_A}")
    print(f"gamma_A = 1 / sqrt(1 - {beta_A}^2) = {gamma_A:.5f}")
    print(f"1/sqrt(2) = {term_1_over_sqrt_2:.5f}")
    print(f"tan(theta_C) = {numerator:.5f} / ({gamma_A:.5f} * ({term_1_over_sqrt_2:.5f} + {beta_A}))")
    print(f"tan(theta_C) = {numerator:.5f} / {denominator:.5f}")
    print(f"tan(theta_C) = {tan_theta_C:.5f}")
    
    print("\nFinally, the angle is calculated by taking the arctangent:")
    print(f"theta_C = arctan({tan_theta_C:.5f})")
    print(f"\nThe angle between the direction of C and that of its parent particle A is {angle_deg:.3f} degrees.")

# Execute the function
calculate_decay_angle()