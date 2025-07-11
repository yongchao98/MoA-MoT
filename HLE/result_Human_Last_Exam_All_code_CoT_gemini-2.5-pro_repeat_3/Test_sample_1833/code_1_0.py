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

    # Convert angle from degrees to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)

    # Calculate sin and cos of the angle
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # For the special case m = M = 1 kg, the formula for tension T is:
    # T = g * sin(θ) * (5 + cos²(θ)) / (1 + cos²(θ))²

    # Calculate the terms in the formula
    cos_theta_sq = cos_theta**2
    numerator = 5 + cos_theta_sq
    denominator = (1 + cos_theta_sq)**2

    # Calculate the final tension
    tension = g * sin_theta * numerator / denominator

    # Output the explanation and the equation with numerical values
    print("The simplified formula for tension (T) when m = M = 1 kg is:")
    print("T = g * sin(θ) * (5 + cos²(θ)) / (1 + cos²(θ))²\n")
    print("Plugging in the values:")
    print(f"g = {g} m/s²")
    print(f"θ = {theta_deg}°\n")
    
    print("The equation with the values is:")
    print(f"T = {g} * sin({theta_deg}°) * (5 + cos²({theta_deg}°)) / (1 + cos²({theta_deg}°))²")
    print(f"T = {g} * {sin_theta:.4f} * (5 + {cos_theta_sq:.4f}) / (1 + {cos_theta_sq:.4f})²")
    print(f"T = {g} * {sin_theta:.4f} * {numerator:.4f} / {denominator:.4f}")

    # Print the final result
    print(f"\nThe calculated tension is: {tension:.2f} N")

calculate_tension()
<<<28.52>>>