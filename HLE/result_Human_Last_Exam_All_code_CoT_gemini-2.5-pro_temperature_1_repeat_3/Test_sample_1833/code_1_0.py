import math

def calculate_tension():
    """
    Calculates the tension in the string for the given physics problem.
    """
    # Given values
    m = 1.0  # kg
    M = 1.0  # kg
    g = 9.8  # m/s^2
    theta_deg = 60.0  # degrees

    # Convert angle to radians for use in math functions
    theta_rad = math.radians(theta_deg)

    # Calculate sin(theta) and cos(theta)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    cos_theta_sq = cos_theta**2

    # The derived formula for tension when m = M is:
    # T = M * g * sin(theta) * (5 + cos(theta)^2) / (1 + cos(theta)^2)^2
    
    # Calculate the tension
    numerator = M * g * sin_theta * (5 + cos_theta_sq)
    denominator = (1 + cos_theta_sq)**2
    tension = numerator / denominator

    # Print the equation with the numerical values
    print("The formula for tension (T) in this specific case (m=M) is:")
    print("T = M * g * sin(θ) * (5 + cos(θ)²) / (1 + cos(θ)²)²\n")
    print("Plugging in the values:")
    # Using format to control the precision for clarity
    print(f"T = {M:.1f} * {g:.1f} * sin({theta_deg:.0f}°) * (5 + cos({theta_deg:.0f}°)²) / (1 + cos({theta_deg:.0f}°)²)²")
    print(f"T = {M:.1f} * {g:.1f} * {sin_theta:.3f} * (5 + {cos_theta_sq:.2f}) / (1 + {cos_theta_sq:.2f})²")
    print(f"T = {M * g * sin_theta:.3f} * {5 + cos_theta_sq:.2f} / ({1 + cos_theta_sq:.2f})²")
    print(f"T = {numerator:.3f} / {denominator:.4f}")
    
    # Print the final result rounded to two decimal places
    print(f"\nThe calculated tension is: {tension:.2f} N")

calculate_tension()