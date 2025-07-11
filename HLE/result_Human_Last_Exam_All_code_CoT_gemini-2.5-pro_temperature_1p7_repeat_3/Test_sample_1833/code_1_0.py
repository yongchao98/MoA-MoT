import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # Given parameters
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    g = 9.8  # acceleration due to gravity in m/s^2
    theta_deg = 60.0  # angle in degrees

    # Convert angle to radians for math functions
    theta_rad = math.radians(theta_deg)

    # Calculate sin and cos of the angle
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # The general formula for tension T is:
    # T = M * m * g * sin(theta) * (3*m + M*(2 + cos(theta)**2)) / (m + M*cos(theta)**2)**2

    # We will print the calculation step-by-step with the given values.
    print("The final equation for tension T is:")
    print("T = (M * m * g * sin(theta)) * (3*m + M*(2 + cos(theta)^2)) / (m + M*cos(theta)^2)^2")
    print("\nPlugging in the values:")
    print(f"M = {M} kg, m = {m} kg, g = {g} m/s^2, theta = {theta_deg} degrees")
    print(f"sin({theta_deg}) = {sin_theta:.4f}")
    print(f"cos({theta_deg}) = {cos_theta:.4f}")
    
    val_m = m
    val_M = M
    val_g = g
    val_sin_theta = sin_theta
    val_cos_theta = cos_theta
    val_cos_theta_sq = val_cos_theta**2

    # Breaking down the calculation for clarity
    numerator_part1 = val_M * val_m * val_g * val_sin_theta
    numerator_part2 = 3 * val_m + val_M * (2 + val_cos_theta_sq)
    denominator = (val_m + val_M * val_cos_theta_sq)**2

    print(f"\nT = ({val_M} * {val_m} * {val_g} * {val_sin_theta:.4f}) * (3*{val_m} + {val_M}*(2 + {val_cos_theta:.4f}^2)) / ({val_m} + {val_M}*{val_cos_theta:.4f}^2)^2")
    print(f"T = ({numerator_part1:.4f}) * (3 + (2 + {val_cos_theta_sq:.4f})) / (1 + {val_cos_theta_sq:.4f})^2")
    print(f"T = {numerator_part1:.4f} * ({numerator_part2:.4f}) / ({denominator:.4f})")

    # Final calculation
    T = (numerator_part1 * numerator_part2) / denominator
    
    print(f"T = {(numerator_part1 * numerator_part2):.4f} / {denominator:.4f}")
    print(f"\nThe calculated tension is {T:.4f} N.")
    print(f"Rounded to two decimal places, the tension is {T:.2f} N.")

# Execute the function
calculate_tension()
<<<28.52>>>