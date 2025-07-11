import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # Given values
    m = 1.0  # mass of the ring in kg
    M = 1.0  # mass of the object in kg
    g = 9.8  # acceleration due to gravity in m/s^2
    theta_deg = 60.0  # angle in degrees

    # Convert angle to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)

    # Calculate sin(theta) and cos(theta)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # The derived formula for tension T
    # T = (M * m * g * sin(theta) * (3*m + M*(2 + cos(theta)^2))) / (m + M*cos(theta)^2)^2

    # Breaking down the calculation
    numerator = M * m * g * sin_theta * (3 * m + M * (2 + cos_theta**2))
    denominator = (m + M * cos_theta**2)**2
    
    tension = numerator / denominator

    # Print the equation with the values substituted
    print("The tension T is calculated using the formula:")
    print("T = (M * m * g * sin(θ) * (3*m + M*(2 + cos(θ)²))) / (m + M*cos(θ)²)²\n")
    print("Substituting the given values:")
    print(f"M = {M} kg")
    print(f"m = {m} kg")
    print(f"g = {g} m/s²")
    print(f"θ = {theta_deg}°\n")
    
    print("The calculation is:")
    print(f"T = ({M} * {m} * {g} * sin({theta_deg}°) * (3*{m} + {M}*(2 + cos({theta_deg}°)^2))) / (({m} + {M}*cos({theta_deg}°)^2)^2)")
    print(f"T = ({M} * {m} * {g} * {sin_theta:.4f} * (3*{m} + {M}*(2 + {cos_theta**2:.4f}))) / (({m} + {M}*{cos_theta**2:.4f})^2)")
    print(f"T = {numerator:.4f} / {denominator:.4f}")
    
    # Print the final result
    print(f"\nThe tension in the string is: {tension:.2f} N")

calculate_tension()