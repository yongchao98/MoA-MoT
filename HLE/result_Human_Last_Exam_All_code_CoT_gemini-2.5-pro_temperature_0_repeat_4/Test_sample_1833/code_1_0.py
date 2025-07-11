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

    # Calculate sin(theta) and cos(theta)^2
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    cos_theta_sq = cos_theta**2

    # The derived formula for tension T is:
    # T = (m*M*g*sin(theta) * (3m + 2M + M*cos(theta)^2)) / (m + M*cos(theta)^2)^2
    
    # Calculate the numerator and denominator parts
    numerator_term = 3 * m + 2 * M + M * cos_theta_sq
    denominator_term = (m + M * cos_theta_sq)**2
    
    numerator = m * M * g * sin_theta * numerator_term
    denominator = denominator_term
    
    # Calculate the final tension
    tension = numerator / denominator

    # Print the final equation with all the numbers substituted
    print("The final equation with the numerical values is:")
    # We use f-string formatting to display the equation and the result.
    # Values are formatted for clarity.
    print(f"T = ({m} * {M} * {g} * {sin_theta:.4f} * (3*{m} + 2*{M} + {M}*{cos_theta_sq:.2f})) / ({m} + {M}*{cos_theta_sq:.2f})^2 = {tension:.2f} N")

calculate_tension()
<<<28.52>>>